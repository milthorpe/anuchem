#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>

#define MAXL (5+1)
#define cart(l) ((l+1)*(l+2)/2)

FILE *fptr; int a,b,cona,conb,n,l; double time,ntime,fwcount;

double min[MAXL][MAXL],Min,
       max[MAXL][MAXL],Max,
       s[MAXL][MAXL],S,
       ss[MAXL][MAXL],SS,
       mean,sd;

long count[MAXL][MAXL],Count;
int F123[MAXL][MAXL],F4[MAXL][MAXL],
    W123[MAXL][MAXL],W4[MAXL][MAXL];

void clear(){
    Min=DBL_MAX; Max=DBL_MIN;
    S=0.; SS=0.;
    Count=0;
    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        min[a][b]=DBL_MAX; max[a][b]=DBL_MIN;
        s[a][b]=0.; ss[a][b]=0.;
        count[a][b]=0;
    }
}

void tally(){
    if (ntime<min[a][b]) min[a][b]=ntime;  if (ntime<Min) Min=ntime;
    if (ntime>max[a][b]) max[a][b]=ntime;  if (ntime>Max) Max=ntime;
    s[a][b]+=ntime; S+=ntime;
    ss[a][b]+=ntime*ntime; SS+=ntime*ntime;
    count[a][b]++; Count++;
}

void printstat(){
    printf(" a  b       n          min         mean         s.d.          max\n-----------------------------------------------------------------\n");
    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        mean=s[a][b]/count[a][b];
        sd=sqrt(ss[a][b]/count[a][b]-mean*mean);
        if (count[a][b]) printf("%2d %2d %7ld %e %e %e %e\n",a,b,count[a][b],min[a][b],mean,sd,max[a][b]);
    }
    mean=S/Count; sd=sqrt(SS/Count-mean*mean);
    printf("Total %7ld %e %e %e %e\n\n",Count,Min,mean,sd,Max);
}

int main(int argc, char *argv[]) {
    if (argc<2) {printf("Please specify output file.\n"); exit(1);}

    // READ FETCH/WRITE COST PARAMETERS
    fptr=fopen("cost.param","r");
    //printf(" a  b  F123    F4  W123    W4\n-------------------------------\n");
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d",&a,&b);
        fscanf(fptr,"%d %d %d %d",&F123[a][b],&F4[a][b],&W123[a][b],&W4[a][b]);
        //printf("%2d %2d %5d %5d %5d %5d\n",a,b,F123[a][b],F4[a][b],W123[a][b],W4[a][b]);
    }
    fclose(fptr);

    // TALLY TIME PER FETCH & WRITE
    if (!(fptr=fopen(argv[1],"r")))  {printf("Can not open '%s'.\n",argv[1]); exit(1);} else printf("Reading '%s'...\n",argv[1]);
    clear();
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d %d %d %d %d %lf",&a,&b,&cona,&conb,&n,&l,&time);
        fwcount=pow(l+1.,2.)*((F123[a][b]+W123[a][b])*cona*conb+F4[a][b]+W4[a][b]);
        ntime=time/fwcount;   
        tally();     
    }
    fclose(fptr);
    printf("Timing statistics per FETCH/WRITE\n");
    printstat();

    // TALLY TIME PER CONTRACTED INTEGRAL
    if (!(fptr=fopen(argv[1],"r")))  {printf("Can not open '%s'.\n",argv[1]); exit(1);} else printf("Reading '%s'...\n",argv[1]);
    clear();
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d %d %d %d %d %lf",&a,&b,&cona,&conb,&n,&l,&time);
        ntime=time/pow(l+1.,2.)/cart(a)/cart(b);
        tally();
    }
    fclose(fptr); 
    printf("Timing statistics per FETCH/WRITE\n");
    printstat();

    // TALLY TIME PER GENCLASS CALL
    if (!(fptr=fopen(argv[1],"r")))  {printf("Can not open '%s'.\n",argv[1]); exit(1);} else printf("Reading '%s'...\n",argv[1]);
    clear();
    while (!feof(fptr)) {
        fscanf(fptr,"%d %d %d %d %d %d %lf",&a,&b,&cona,&conb,&n,&l,&time);
        ntime=time;
        tally();
    }
    fclose(fptr);
    printf("Timing statistics per FETCH/WRITE\n");
    printstat();

}
