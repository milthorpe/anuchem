#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define MAXL (5+1)
#define cart(l) ((l+1)*(l+2)/2)
#define read fscanf(fptr,"%d %d %d %d %d %d %lf",&a,&b,&cona,&conb,&n,&l,&time);

FILE *fptr;  char fname[256];
double min[MAXL][MAXL],Min,
       max[MAXL][MAXL],Max,
       s[MAXL][MAXL],S,
       ss[MAXL][MAXL],SS,
       time,mean,sd;
long long count[MAXL][MAXL],Count;
int a,b,cona,conb,n,l,w,
    F123[MAXL][MAXL],F4[MAXL][MAXL],W123[MAXL][MAXL],W4[MAXL][MAXL];

void init(){
    if (!(fptr=fopen(fname,"r")))  {printf("Can not open '%s'.\n",fname); exit(2);} else printf("Reading '%s'...\n",fname);
    Min=DBL_MAX; Max=DBL_MIN; S=0.; SS=0.; Count=0;
    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        min[a][b]=DBL_MAX; max[a][b]=DBL_MIN;
        s[a][b]=0.; ss[a][b]=0.; count[a][b]=0;
    }
}

void tally(){
    double ntime=time/w,sntime=w*ntime*ntime; 
    if (ntime<min[a][b]) min[a][b]=ntime;  if (ntime<Min) Min=ntime;
    if (ntime>max[a][b]) max[a][b]=ntime;  if (ntime>Max) Max=ntime;
    s[a][b]+=time; S+=time;
    ss[a][b]+=sntime; SS+=sntime;
    count[a][b]+=w; Count+=w;
}

void printstat(char *title){
    fclose(fptr); if (!Count) {printf("No data in the file.\n"); exit(3);}
    printf("%s\n",title);
    printf(" a  b                 n          min         mean         s.d.          max\n----------------------------------------------------------------------------\n");
    for (a=0; a<MAXL; a++) for (b=0; b<=a; b++) {
        mean=s[a][b]/count[a][b]; sd=sqrt(ss[a][b]/count[a][b]-mean*mean);
        if (count[a][b]) printf("%2d %2d %17lld %e %e %e %e\n",a,b,count[a][b],min[a][b],mean,sd,max[a][b]);
    }
    mean=S/Count; sd=sqrt(SS/Count-mean*mean);
    printf("Total %17lld %e %e %e %e\n\n",Count,Min,mean,sd,Max);
}

int main(int argc, char *argv[]) {
    if (argc<2) {printf("Please specify output file.\n"); exit(1);} else strcpy(fname,argv[1]);

    // TALLY X PER CONTRACTED INTEGRAL
    init();
    while (!feof(fptr)) {
        read;
        w=(l+1)*(l+1)*cart(a)*cart(b);
        tally();
    }
    printstat("X statistics per CONTRACTED INTEGRAL");

    // TALLY X PER INTPACK CALL
    init(); w=1;
    while (!feof(fptr)) {
        read;
        tally();
    }
    printstat("X statistics per INTPACK CALL");

    printf("Total X = %e\n",S);
}
