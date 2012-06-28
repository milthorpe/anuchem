#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#define MAXL (5+1)
#define cart(l) ((l+1)*(l+2)/2)

double min[MAXL][MAXL],
       max[MAXL][MAXL],
       //tot[MAXL][MAXL],
       s[MAXL][MAXL],
       ss[MAXL][MAXL];

int count[MAXL][MAXL],
    F123[MAXL][MAXL],F4[MAXL][MAXL],
    W123[MAXL][MAXL],W4[MAXL][MAXL];

int main(int argc, char *argv[]) {
    FILE *fptr; int a,b,cona,conb,n,l; double time,ntime;

    if (argc<2) {printf("Please specify output file.\n"); exit(1);}

    fptr=fopen("cost.param","r");
    printf(" a  b  F123    F4  W123    W4\n-------------------------------\n");
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d",&a,&b);
        fscanf(fptr,"%d %d %d %d",&F123[a][b],&F4[a][b],&W123[a][b],&W4[a][b]);
        printf("%2d %2d %5d %5d %5d %5d\n",a,b,F123[a][b],F4[a][b],W123[a][b],W4[a][b]);
    }
    fclose(fptr);

    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        min[a][b]=DBL_MAX;
        max[a][b]=DBL_MIN;
        s[a][b]=0.;
        ss[a][b]=0.;
        count[a][b]=0;
    }

    fptr=fopen(argv[1],"r");
    if (fptr==NULL)  {printf("Can not open '%s'.\n",argv[1]); exit(1);}
    else printf("Reading '%s'...\n",argv[1]);
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d %d %d %d %d %lf",&a,&b,&cona,&conb,&n,&l,&time);
        //printf("%d %d %d %d %d %d %e",a,b,cona,conb,n,l,time);
        ntime=time/pow(l+1.,2.)/cart(a)/cart(b);
        if (ntime<min[a][b]) min[a][b]=ntime;
        if (ntime>max[a][b]) max[a][b]=ntime;
        s[a][b]+=ntime;
        ss[a][b]+=ntime*ntime;
        count[a][b]++;
    }
    fclose(fptr);
    //printf("%d %d %d %d %d %d %e",a,b,cona,conb,n,l,time);
    printf("Sucessfully read '%s'.\n",argv[1]);

    printf(" a  b      n          min         mean         s.d.          max\n-----------------------------------------------------------------\n");
    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        double mean=s[a][b]/count[a][b];
        double sd=sqrt(ss[a][b]/count[a][b]-mean*mean);
        if (count[a][b]) printf("%2d %2d %6d %e %e %e %e\n",a,b,count[a][b],min[a][b],mean,sd,max[a][b]);
    }

}
