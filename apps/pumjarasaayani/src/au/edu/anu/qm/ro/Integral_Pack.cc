#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bessel4.h"
#include "Integral_Pack.h"
#include "cblas.h"

#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

int delta[3][3]={{1,0,0},{0,1,0},{0,0,1}};
double JpY00[11]={0.28209479177387814, 0.09403159725795937, 0.018806319451591877,
       		0.0026866170645131254, 0.000298513007168125,0.000027137546106193183, 
                2.087503546630245e-6,1.39166903108683e-7, 8.186288418157823e-9,
                4.308572851662012e-10, 2.0517013579342914e-11}; 
// 0-10 accomodate (hh|phi) - need more for higer angular momentum - see RO#7 eqn 6b 

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
        printf("Integral_Pack.cc N=%d L=%d\n",N,L);
    }

    Integral_Pack::~Integral_Pack() {
        free(lambda);
        free(q);
        for (int i=1; i<=MAX_BRA_L; i++) for (int j=1; j<=i; j++) {
            free(HRRMAP[i][j]);
        }
    }

    void Integral_Pack::initialize() {
        int x,y,z,l,m=0;
        totalBraL[0]=0;
        totalBraL[1]=1;
        for (l=0; l<=MAX_BRA_L; l++) {
            noOfBra[l]=((l+1)*(l+2))/2;
            if (l>0) totalBraL[l+1]=totalBraL[l]+noOfBra[l];
              for (x=0; x<=l; x++) for (y=0; y<=l-x; y++) { // must be consistent with PowerList.x10

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

             	int aInt=totalBraL[i]+a,bInt=totalBraL[j]+b;
             	int ax=inverseMap3[aInt].x, ay=inverseMap3[aInt].y, az=inverseMap3[aInt].z,
             		bx=inverseMap3[bInt].x,	by=inverseMap3[bInt].y,	bz=inverseMap3[bInt].z;

             	int increment;
             	if (bx) increment=0; else if (by) increment=1; else /*if (bz)*/ increment=2;
             	int bm1 = map3[bx-delta[0][increment]][by-delta[1][increment]][bz-delta[2][increment]]-totalBraL[j-1];
             	int ap1 = map3[ax+delta[0][increment]][ay+delta[1][increment]][az+delta[2][increment]]-totalBraL[i+1];

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

    int Integral_Pack::Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, double* temp, int n, int L){
        int bra,K=(L+1)*(L+1),p,e,i,ii,j,jj,k,l,m; 
        bool swap,lzero=(lambda[n]==0.); // should be replace by lambda[n]>1e-y 
        double (*V1)[K]=(double (*)[K])malloc(totalBraL[a+b+1]*K*sizeof(double));
        double (*V2)[K]=(double (*)[K])malloc(totalBraL[a+b+1]*K*sizeof(double));
        if (V1==NULL || V2==NULL) {printf("Integral_Pack.cc V1/V2 allocation failed size=%d*sizeof(double)\n",totalBraL[a+b+1]*K); exit(1);}
        double (*HRR[MAX_BRA_L+1][MAX_BRA_L+1])[K];
        for (i=a; i<=a+b; i++) {
            HRR[i][0] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]));
            if (HRR[i][0]==NULL) {printf("Integral_Pack.cc HRR[%d][0] allocation failed sized=%d*sizeof(double)\n",i,K*noOfBra[i]); exit(1);}
            memset(HRR[i][0],0,sizeof(double)*K*noOfBra[i]);
        }
        double J[L+a+b+1], Y[(L+1)*(L+1)], rAB2=sqr(A[0]-B[0])+sqr(A[1]-B[1])+sqr(A[2]-B[2]);
        double onelambda;
        if (!lzero) onelambda=-1.0/lambda[n];

        // printf("A %e %e %e B %e %e %e\n",A[0],A[1],A[2],B[0],B[1],B[2]); 
        for (ii=0; ii<dconA; ii++) for (jj=0; jj<dconB; jj++) {
            double zeta=zetaA[ii]+zetaB[jj];
            double P[3]={(zetaA[ii]*A[0]+zetaB[jj]*B[0])/zeta,(zetaA[ii]*A[1]+zetaB[jj]*B[1])/zeta,(zetaA[ii]*A[2]+zetaB[jj]*B[2])/zeta};
            double gAB=exp(-zetaA[ii]*zetaB[jj]/zeta*rAB2)*pow(PI/zeta,1.5)*conA[ii]*conB[jj];
            double one2zeta=.5/zeta;

            double r=sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
            bool rzero=(r==0.); // should be replace by r>1e-y 

            //printf("ii=%d jj=%d zetaA=%e zetaB=%e zeta=%e conA=%e conB=%e rAB2=%e gAB=%e\n",ii,jj,zetaA[ii],zetaB[jj],zeta,conA[ii],conB[jj],rAB2,gAB);

            if (!rzero) { 
                double phi=atan2(P[1],P[0]),X=P[2]/r;
                GenY(Y,X,phi,L);
                GenJ(J,r*lambda[n],L+a+b);
            }

            swap = false;
            for (p=a+b; p>=0; p--) {
                swap = !swap;
                double (* Va)[K] = swap ? V2 : V1;
                double (* Vb)[K] = swap ? V1 : V2;
                memset(Va,0,sizeof(double)*K*totalBraL[a+b+1]);

                // lambda[n]=0.0 is taken care of by separately. (3-term RR & trivial initial conditions) only p=0 contributes! RO#7 eqn11
                if (lzero && p==0) {
                    Va[0][0]=JpY00[0]*q[0]*gAB;
                    for (e=1; e<a+b+1; e++) for (i=0; i<noOfBra[e]; i++) {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])]; // Be careful
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        Va[aplusIndex][0]=(P[j]-A[j])*Va[aIndex][0];
                        if (aj>0) Va[aplusIndex][0] += aj*one2zeta*Va[aminusIndex][0];
                    }
                }

                // Fill e=0
                if (!rzero && !lzero) {
                    double nfactor=q[n]*gAB*exp(-.25*sqr(lambda[n])/zeta)*pow(-.5*lambda[n]/zeta/r,p);
                    for (l=0; l<=L; l++) for (m=-l; m<=l; m++) 
                            Va[0][lm2k(l,m)] = nfactor*J[l+p]*Y[lm2k(l,m)]; // eqn (23) //6a
                }
                else if (rzero && !lzero) 
                    Va[0][0] = q[n]*gAB*exp(-.25*sqr(lambda[n])/zeta)*pow(-.5*sqr(lambda[n])/zeta,p)*JpY00[p]; // l=m=0 only //6b   

                // Fill higher e
            	if (!lzero) {
                    for (e=1; e<a+b+1-p; e++) for (i=0; i<noOfBra[e]; i++)  {
                        int aplusIndex = totalBraL[e]+i;
                        int j=buildMap[aplusIndex]; // printf("j=%d\n",j);
                        int x=inverseMap3[aplusIndex].x,y=inverseMap3[aplusIndex].y,z=inverseMap3[aplusIndex].z;
                        int aIndex = map3[x-delta[0][j]][y-delta[1][j]][z-delta[2][j]];
                        int aminusIndex = map3[abs(x-2*delta[0][j])][abs(y-2*delta[1][j])][abs(z-2*delta[2][j])];
                        int aj = delta[0][j]*(x-1) + delta[1][j]*(y-1) + delta[2][j]*(z-1);
                        double paj = P[j]-A[j];

                        for (l=0; l<=L; l++) for (m=-l; m<=l; m++) {
                   	        int lm=lm2k(l,m),
                                k=lm,
                                kxplus=lm2k(l-1,m+1),
                                kxminus=lm2k(l-1,m-1),
                                kyplus=lm2k(l-1,-(m+1)),
                                kyminus=lm2k(l-1,-(m-1)),
                                kzero=lm2k(l-1,m);
                            if (aIndex>=totalBraL[a+b+1] || aIndex<0 || k<0 || k>K) printf("aIndex=%d k=%d\n",aIndex,k);
                            double vapk = paj*Va[aIndex][k]+P[j]*Vb[aIndex][k]; 
                            if (aj>0) vapk += aj*one2zeta*(Va[aminusIndex][k]+Vb[aminusIndex][k]);
                            //printf("[%d %d %d | %2d %2d %2d] = %f ainx=%d\n",inverseMap3[aplusIndex].x,inverseMap3[aplusIndex].y,inverseMap3[aplusIndex].z,
                            //		n,l,m,vapk,aIndex);

                            if (l!=0) {
                                // no extra term for l=0 
                                if (j==2) vapk += onelambda*cz[lm]*Vb[aIndex][kzero];
                                else if (j==1) vapk += onelambda*(cyminus[lm]*Vb[aIndex][kyminus]+cyplus[lm]*Vb[aIndex][kyplus]);
                                else /*if (j==0)*/ vapk += onelambda*(cxminus[lm]*Vb[aIndex][kxminus]+cxplus[lm]*Vb[aIndex][kxplus]);
                                        // It's better to get the if statement out of the loop
                                        //printf("[%d %d %d | %2d %2d %2d] = %.15e j=%d cy+ =%e aj=%d\n",inverseMap3[aplusIndex].x,inverseMap3[aplusIndex].y,inverseMap3[aplusIndex].z,
                                       //		n,l,m,vapk,j,cxminus[lm]*Vb[aIndex][kxminus], aj );
                            }
                            if (aplusIndex>=totalBraL[a+b+1] || aplusIndex<0 || k<0 || k>K) printf("aplusIndex=%d k=%d\n",aplusIndex,k);
                            Va[aplusIndex][k] = vapk; 
                        }
                    }
                }

            }
            double (* Va)[K] = swap ? V2 : V1;
            for (i=a; i<=a+b; i++) for (bra=0; bra<noOfBra[i]; bra++) 
                cblas_daxpy(K, 1.0, Va[bra+totalBraL[i]], 1, HRR[i][0][bra], 1);
        }

        free(V1);free(V2);

        double dd[3]={A[0]-B[0],A[1]-B[1],A[2]-B[2]};
   
        for (j=1; j<=b; j++) for (i=a; i<=a+b-j; i++)  {
            HRR[i][j] = (double (*)[K])(malloc(sizeof(double)*K*noOfBra[i]*noOfBra[j]));
            if (HRR[i][j]==NULL) {printf("Integral_Pack.cc HRR[%d][%d] size=%d*sizeof(double)\n",i,j,K*noOfBra[i]*noOfBra[j]); exit(1);}
            for (ii=0; ii<noOfBra[i]; ii++) for (jj=0; jj<noOfBra[j]; jj++) {
            	int lindex = ii*noOfBra[j] + jj;
            	int rindex1=HRRMAP[i][j][lindex].x;
            	int rindex2=HRRMAP[i][j][lindex].y;

                double* lhs = HRR[i][j][lindex];
                double* rhs1 = HRR[i+1][j-1][rindex1];
                double* rhs2 = HRR[i][j-1][rindex2];

            	double factor = dd[HRRMAP[i][j][lindex].z];

                cblas_dcopy(K, rhs1, 1, lhs, 1);
                cblas_daxpy(K, factor, rhs2, 1, lhs, 1);
        	}
            //if (HRR[i][j-1]==NULL) {printf("Integral_Pack.cc ln337\n"); exit(1);}
            free(HRR[i][j-1]);
            if (i==a+b-j) {
                //if (HRR[i+1][j-1]==NULL) {printf("Integral_Pack.cc ln340\n"); exit(1);}
                free(HRR[i+1][j-1]);
            }
        }

        int totalInt = noOfBra[a]*noOfBra[b]*K;
        memcpy(temp, HRR[a][b], totalInt*sizeof(double));
        //if (HRR[a][b]==NULL) {printf("Integral_Pack.cc ln362\n"); exit(1);}
        free(HRR[a][b]);
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
