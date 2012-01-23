#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bessel4.h"
#include "Integral_Pack.h"

#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

int delta[3][3]={{1,0,0},{0,1,0},{0,0,1}};

namespace au {
    namespace edu {
        namespace anu {
            namespace qm {
                namespace ro {

    Integral_Pack* Integral_Pack::_make(int N, int L) {
      return new Integral_Pack(N,L);
    }

    Integral_Pack::Integral_Pack(int N,int L) {
        this->N = N; this->L=L;
        initialize();
        initializeCoulomb(N);
    }

    void Integral_Pack::initialize() {
        int x,y,z,l,m=0;
        totalBraL[0]=1;
        for (l=0; l<=MAX_BRA_L; l++) {
            noOfBra[l]=((l+1)*(l+2))/2;
            if (l>0) totalBraL[l]=totalBraL[l-1]+noOfBra[l];
            for (x=l; x>=0; x--) for (y=l-x; y>=0; y--) {

                map3[x][y][z=l-x-y] = m;
                inverseMap3[m].x = x;
                inverseMap3[m].y = y;
                inverseMap3[m].z = z;

                // LHG 2012 Table II & boolean before (36)
                if (z == 1) buildMap[m]=2; //z
                else if (y == 1) buildMap[m]=1; //y
                else if (x == 1) buildMap[m]=0; //x
                else if (z >= 2) buildMap[m]=2; //z
                else if (y >= 2) buildMap[m]=1; //y
                else buildMap[m]=0; //x

                //printf("%d : %d %d %d j=%d\n", m, x,y,z, buildMap[m]);
                m++;
            }
        }

        // HRR
        int i,j,a,b;
        for (i=1; i<=MAX_BRA_L; i++) for (j=1; j<=i; j++) {
            HRRMAP[i][j] = (Point *) malloc(sizeof(Point)*noOfBra[i]*noOfBra[j]);
            for (a=0; a<noOfBra[i]; a++) for (b=0; b<noOfBra[j]; b++) {

             	int aInt=totalBraL[i-1]+a,bInt=totalBraL[j-1]+b;
             	int ax=inverseMap3[aInt].x, ay=inverseMap3[aInt].y, az=inverseMap3[aInt].z,
             		bx=inverseMap3[bInt].x,	by=inverseMap3[bInt].y,	bz=inverseMap3[bInt].z;

             	int increment;
             	if (bx) increment=0; else if (by) increment=1; else /*if (bz)*/ increment=2;
             	int bm1 = map3[bx-delta[0][increment]][by-delta[1][increment]][bz-delta[2][increment]]-(j>1?totalBraL[j-2]:0);
             	int ap1 = map3[ax+delta[0][increment]][ay+delta[1][increment]][az+delta[2][increment]]-totalBraL[i];

             	int leftindex=noOfBra[j]*a+b;
             	int rightindexA=noOfBra[j-1]*ap1+bm1;
             	int rightindexB=noOfBra[j-1]*a+bm1;
                HRRMAP[i][j][leftindex].x=rightindexA;
                HRRMAP[i][j][leftindex].y=rightindexB;
                HRRMAP[i][j][leftindex].z=increment;
            }
        }

        //printf("%d\n",MAX_TOTAL_BRA_L);
        // LMR 2012 manuscript in preparation - similar to LHG 2012 eqn (27)
        for (l=0; l<=MAX_KET_L; l++) for (m=-l; m<=l; m++) {
            double cminus=.5*sqrt((l-m-1.0)*(l-m)*(2.0*l+1.0)/(2.0*l-1.0)),
                cplus=.5*sqrt((l+m-1.0)*(l+m)*(2.0*l+1.0)/(2.0*l-1.0));
            int index=lm2k(l,m);

            if (m<=-2) cxplus[index]=-cminus;
                else if (m>=1) cxplus[index]=cminus;
                else if (m==-1) cxplus[index]=0;
                else /*if (m==0)*/ cxplus[index]=SQRT2*cminus;

            if (m<=-1) cxminus[index]=cplus;
                else if (m>=2) cxminus[index]=-cplus;
                else if (m==0) cxminus[index]=0;
                else /*if (m==1)*/ cxminus[index]=-SQRT2*cplus;

            if (m<=-2) cyplus[index] = -cminus;
                else if (m>=1) cyplus[index] = cminus;
                else if (m==-1) cyplus[index] = -SQRT2*cminus;
                else /*if (m==0)*/ cyplus[index] = SQRT2*cminus;

            if (m<=-1) cyminus[index]=-cplus;
                else if (m>=2) cyminus[index]=cplus;
                else /*if (m==0)*/ cyminus[index]=0;
                //else if (m==1) cyminus[index]=0;

            cz[index] = sqrt((l*l-m*m)*(2.0*l+1.0)/(2.0*l-1.0));
            //printf("%2d %2d : %15.6e %15.6e %15.6e %15.6e %15.6e : %15.6e %15.6e %d\n",l,m,cxplus[index],cxminus[index],cyplus[index],cyminus[index],cz[index],cplus,cminus,index);
        }
    }

    int Integral_Pack::GenJ(double *B, double x, int L) {
        int l,nb=L+1,ncalc;

        if (x<1e-20) { // To guard division by zero later in "fac" - 1e-20 should be justified by checking with Mathematica.
            for (l=0; l<=L; l++) B[l]=0.;
            B[0]=1.;
            return 0;
        }

        J_bessel_HalfInt(x, nb, B, &ncalc);
        //if (ncalc!=L+1) printf("bessel %d %f\n",ncalc,x); // Check what happens? Zero out the rest?
        double fac = sqrt(PI*.5/x); 
        for (l=0; l<ncalc; l++) B[l]*=fac;
        for (; l<=L; l++) B[l]=0;
    }

    int Integral_Pack::GenY(double *Y, double X, double phi, int L) {
        // Plm calculation according to (6.7.9)-(6.7.10) NR 3rd ed
        int l,m;
        double Plm[L+1][L+1],sintheta = sqrt(1-X*X);

        Plm[0][0] = 0.5/sqrt(PI);
        for (l=1;l<=L; l++)
            Plm[l][l]=-sqrt(1.0+0.5/l)*sintheta*Plm[l-1][l-1];

        for (m=0; m<L; m++) Plm[m+1][m]=X*sqrt(2.0*m+3.0)*Plm[m][m];

        for (l=2; l<=L; l++) {
            double ls=l*l, lm1s = (l-1)*(l-1);
            for (m=0; m <= l-2; m++) {
                double ms=m*m;
                Plm[l][m] = sqrt((4.0*ls-1.0)/(ls-ms))*(X*Plm[l-1][m]-sqrt((lm1s-ms)/(4.0*lm1s-1.0))*Plm[l-2][m]);
            }
        }

        // Real Ylm
        for (l=0; l<=L; l++) {
            Y[lm2k(l,0)]  = Plm[l][0];
            for (m=1; m<=l; m++) {
                Y[lm2k(l,m)]  = Plm[l][m]*SQRT2*cos(m*phi);
                Y[lm2k(l,-m)] = Plm[l][m]*SQRT2*sin(m*phi);
            }
        }
    }

    int Integral_Pack::Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp){
        int bra,K = (N+1)*(L+1)*(L+1),p,e,i,ii,j,jj,k,n,l,m,initialn=lambda[0]==0.?1:0; // more regorious check / cut-off required
        double V1[totalBraL[a+b]][K];
        double V2[totalBraL[a+b]][K];
        bool swap = false;
        double V[totalBraL[a+b]][K];
        memset(V,0,totalBraL[a+b]*K*sizeof(double));

        double Y[(L+1)*(L+1)],J[(L+a+b+1)*(N+1)],JpY00[9]={0.28209479177387814, 0.09403159725795937, 0.018806319451591877,
        		0.0026866170645131254, 0.000298513007168125,0.000027137546106193183, 2.087503546630245e-6,1.39166903108683e-7, 8.186288418157823e-9};

        for (ii=0; ii<dconA; ii++) for (jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta, (zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta, (zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double rAB2 = sqr(A[0]-B[0])+sqr(A[1]-B[1])+sqr(A[2]-B[2]);
            double gAB=exp(-zetaA[ii]*zetaB[jj]/zeta*rAB2)*pow(PI/zeta,1.5)*conA[ii]*conB[jj];
            double one2zeta = .5/zeta;
            printf("ii=%d jj=%d zetaA=%e zetaB=%e zeta=%e conA=%e conB=%e rAB2=%e gAB=%e\n",ii,jj,zetaA[ii],zetaB[jj],zeta,conA[ii],conB[jj],rAB2,gAB);
            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            double X=P[2]/r;
            double phi=atan2(P[1],P[0]);

            if (r!=0.0) { // should be replace by r>1e-y 
                GenY(Y,X,phi,L);
                for (n=0; n<=N; n++) GenJ(&J[(L+a+b+1)*n],r*lambda[n],L+a+b);
            }

            for (p=a+b; p>=0; p--) {
                swap = !swap;
                double (*Va)[totalBraL[a+b]][K] = swap ? &V2 : &V1;
                double (*Vb)[totalBraL[a+b]][K] = swap ? &V1 : &V2;

                // lambda[0]=0 is taken care of by separately. (3-term RR & trivial initial conditions)
        	    // only p=0 contributes!
                if (initialn==1 && p==0) {
                    memset(Va,0,sizeof(double)*K*totalBraL[a+b]);
                    (*Va)[0][0]=JpY00[0]*q[0]*gAB;

                    int nOffset=0;
                    for (e=0; e<a+b; e++) for (i=0; i<noOfBra[e+1]; i++) {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])]; // Be careful
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        (*Va)[aplusIndex][nOffset]=(P[j]-A[j])*(*Va)[aIndex][nOffset];
                        if (aj>0) (*Va)[aplusIndex][nOffset] += aj*one2zeta*(*Va)[aminusIndex][nOffset];
                    }
                }

                for (n=initialn; n<=N; n++) {
                    if (r!=0.0) {
                        double nfactor=q[n]*gAB*exp(-.25*sqr(lambda[n])/zeta)*pow(-.5*lambda[n]/zeta/r,p);
                        for (l=0; l<=L; l++) for (m=-l; m<=l; m++) {
                            (*Va)[0][n*(L+1)*(L+1)+lm2k(l,m)] = nfactor*J[(L+a+b+1)*n+l+p]*Y[lm2k(l,m)]; // eqn (23)
                            //printf("[0 0 0 | %2d %2d %2d] ^%d = %15.8e J= %15.8e Y=%15.8e\n",n,l,m,p, (*Va)[0][n*(L+1)*(L+1)+lm2k(l,m)],J[(L+a+b+1)*n+l+p],Y[lm2k(l,m)]);
                        }
                    }
                    else {
                        double nfactor=q[n]*gAB*exp(-.25*sqr(lambda[n])/zeta)*pow(-.5*sqr(lambda[n])/zeta,p);
                        memset(Va,0,sizeof(double)*K*totalBraL[a+b]);
                        (*Va)[0][n*(L+1)*(L+1)] = nfactor*JpY00[p]; // l=m=0 only
                    }
                }


        	    for (n=initialn; n<=N; n++) {
                    int nOffset=n*(L+1)*(L+1);
                    double onelambda=-1.0/lambda[n];

                    for (e=0; e<a+b-p; e++) for (i=0; i<noOfBra[e+1]; i++)  {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])];
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);

                        for (l=0; l<=L; l++) for (m=-l; m<=l; m++) {
                   	        int lm=lm2k(l,m),
                                k=nOffset+lm,
                                kxplus=nOffset+lm2k(l-1,m+1),
                                kxminus=nOffset+lm2k(l-1,m-1),
                                kyplus=nOffset+lm2k(l-1,-(m+1)),
                                kyminus=nOffset+lm2k(l-1,-(m-1)),
                                kzero=nOffset+lm2k(l-1,m);

                            (*Va)[aplusIndex][k]=(P[j]-A[j])*(*Va)[aIndex][k]+P[j]*(*Vb)[aIndex][k];
                            if (aj>0) (*Va)[aplusIndex][k] += aj*one2zeta*((*Va)[aminusIndex][k]+(*Vb)[aminusIndex][k]);

                                    //printf("[%d %d %d | %2d %2d %2d] = %f ainx=%d\n",inverseMap3[aplusIndex].x,inverseMap3[aplusIndex].y,inverseMap3[aplusIndex].z,
                                    //		n,l,m,(*Va)[aplusIndex][k],aIndex);

                            if (l==0) continue;
                                    // get if out of loop over l and m
                            if (j==2) (*Va)[aplusIndex][k] += onelambda*cz[lm]*(*Vb)[aIndex][kzero];
                            else if (j==1) (*Va)[aplusIndex][k] += onelambda*(cyminus[lm]*(*Vb)[aIndex][kyminus]+cyplus[lm]*(*Vb)[aIndex][kyplus]);
                            else /*if (j==0)*/ (*Va)[aplusIndex][k] += onelambda*(cxminus[lm]*(*Vb)[aIndex][kxminus]+cxplus[lm]*(*Vb)[aIndex][kxplus]);
                                    // It's better to get the if statement out of the loop
                                    //printf("[%d %d %d | %2d %2d %2d] = %.15e j=%d cy+ =%e aj=%d\n",inverseMap3[aplusIndex].x,inverseMap3[aplusIndex].y,inverseMap3[aplusIndex].z,
                                   //		n,l,m,(*Va)[aplusIndex][k],j,cxminus[lm]*(*Vb)[aIndex][kxminus], aj );

                        }
                    }
                }

            }
            double (*Va)[totalBraL[a+b]][K] = swap ? &V2 : &V1;
            for (bra=0; bra<totalBraL[a+b]; bra++) for(k=0; k<K; k++) V[bra][k] += (*Va)[bra][k];
        }

        //printf("[ x y z | n l m ] (a=%d b=%d nof=%d) \n",a,b,noOfBra[a+b]);
        bra=-1;
        for (p=0; p<=a+b; p++) for (i=0; i<noOfBra[p]; i++) {
            bra++;
            for (n=0; n<=N; n++) for (l=0; l<=L; l++) for (m=-l; m<=l; m++)
                //printf("[%d %d %d | %2d %2d %2d] = %25.15e\n",inverseMap3[bra].x,inverseMap3[bra].y,inverseMap3[bra].z,n,l,m,(*Va)[bra][n*(L+1)*(L+1)+lm2k(l,m)]);
                printf("%.15e\n",V[bra][n*(L+1)*(L+1)+lm2k(l,m)]);
        }
        // HRR
        double dd[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};
        double (*HRR[MAX_BRA_L+1][MAX_BRA_L+1])[K];
        double (*Va)[totalBraL[a+b]][K] = swap ? &V2 : &V1;
        for (i=a; i<=a+b; i++) {
            HRR[i][0] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]));
            for (bra=0; bra<noOfBra[i]; bra++ ) for (k=0; k<K; k++) {
                HRR[i][0][bra][k] = (*Va)[bra+(i>0?totalBraL[i-1]:0)][k];
            }
        }

        for (i=a; i<=a+b; i++) for (j=1; j<=i-a; j++) {
            HRR[i-j][j] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i-j]*noOfBra[j]));
            for (ii=0; ii<noOfBra[i-j]; ii++) for (jj=0; jj<noOfBra[j]; jj++) {
            	int lindex = ii*noOfBra[j] + jj;
            	int rindex1=HRRMAP[i-j][j][lindex].x;
            	int rindex2=HRRMAP[i-j][j][lindex].y;
            	double factor = dd[HRRMAP[i-j][j][lindex].z];
            	for (k=0; k<K; k++)
                    HRR[i-j][j][lindex][k] = HRR[i-j+1][j-1][rindex1][k]+factor*HRR[i-j][j-1][rindex2][k];
        	}
        }
        int ind=0;
        for (ii=0; ii<noOfBra[a]; ii++) for (jj=0; jj<noOfBra[b]; jj++){
        	int lindex = ii*noOfBra[b] + jj;
        	int indexa = ii + (a>0? totalBraL[a-1]:0);
        	int indexb = jj + (b>0? totalBraL[b-1]:0);
            for (n=0; n<=N; n++) for (l=0; l<=L; l++) for (m=-l; m<=l; m++) {
                printf("[%d,%d,%d %d,%d,%d | %2d %2d %2d] = %25.15e\n",inverseMap3[indexa].x,inverseMap3[indexa].y,inverseMap3[indexa].z,
                		inverseMap3[indexb].x,inverseMap3[indexb].y,inverseMap3[indexb].z,n,l,m,HRR[a][b][lindex][n*(L+1)*(L+1)+lm2k(l,m)]);
                temp[ind++]=HRR[a][b][lindex][n*(L+1)*(L+1)+lm2k(l,m)];
            }
        }

    }

    int Integral_Pack::initializeCoulomb(int N){
        lambda = (double *) malloc(sizeof(double)*(N+1));
        q = (double *) malloc(sizeof(double)*(N+1));
        int i;
        for (i=0; i<=N; i++) {
            lambda[i]=i;
            q[i]=2.*sqrt(2.);
        }
        q[0]=2.;
    }

// end Integral_Pack

                }
            }
        }
    }
}

/*
int main() {

    // Test Ylm
    //int l,m; double Y[MAX_KET_LM],X=0.5941074799900165,phi=0.8002747424752843;
    //GenY(Y,X,phi,MAX_KET_L);
    //for (l=0; l<=MAX_KET_L; l++) {
    //	for (m=-l; m<=l; m++) printf("%8.5e  ",Y[lm2k(l,m)]);
    //	printf("\n");
    //}

    // Test jl
    //int l; double x=.567,J[MAX_KET_L+1];
    //GenJ(J,x,MAX_KET_L);
    //for (l=0; l<=MAX_KET_L; l++)
    //	printf("%8.5e  ",J[l]);
    //printf("\n");

    // Test Genclass
    int N=10,L=10,dconA=6,dconB=3;
    double A[3],B[3];
    double conA[dconA],conB[dconB],zetaA[dconA],zetaB[dconB];

    au::edu::anu::qm::ro::Integral_Pack* ip = new au::edu::anu::qm::ro::Integral_Pack(N);

    A[0]=.1; B[0]=.4;
    A[1]=.2; B[1]=.3;
    A[2]=.3; B[2]=.2;

    zetaA[0]=1.;  conA[0]=6.;
    zetaA[1]=2.;  conA[1]=5.;
    zetaA[2]=3.;  conA[2]=4.;
    zetaA[3]=4.;  conA[3]=3.;
    zetaA[4]=5.;  conA[4]=2.;
    zetaA[5]=6.;  conA[5]=1.;

    zetaB[0]=1.;  conB[0]=1.5;
    zetaB[1]=2.;  conB[1]=2.5;
    zetaB[2]=4.5; conB[2]=3.5;

    // int Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB,
    // double *conA, double *conB, int dconA, int dconB,
    // double *lambda, double *q, int N, int L)
    ip->Genclass(2,0,A,B,zetaA,zetaB,conA,conB,dconA,dconB);

}
*/
