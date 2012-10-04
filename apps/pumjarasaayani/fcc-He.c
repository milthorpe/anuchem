#include<stdio.h>
#include<math.h>
#include<string.h>
#define NND 4.2 // nearest neighbor distance (A) - JCP 129 134107 (2008) pp 9 Tab II

int main (int argc, char *argv[]){
    double x,y,z;
    int i,j,k,n=1;
    if (argc>1) n=atoi(argv[1]);
    while (n) {
        // first layer
        printf("%d\nfcc\n",n*n);
        for (j=0; j<n; j++) // vertical (y-axis)
            for (i=0; i<n-j; i++) { // horizontal (x-axis)
                 x=j*NND*.5+i*NND;
                 y=j*NND*sqrt(3.)/2.;
                 z=0.;
                 printf("He\t%10.3f\t%10.3f\t%10.3f\n",x,y,z);
            }
        // second layer
        n--;
        for (j=0; j<n; j++) // vertical (y-axis)
            for (i=0; i<n-j; i++) { // horizontal (x-axis)
                 x=.5*NND+j*NND*.5+i*NND;
                 y=sqrt(3.)/4.*NND+j*NND*sqrt(3.)/2.;
                 z=NND*sqrt(3.)/2.;
                 printf("He\t%10.3f\t%10.3f\t%10.3f\n",x,y,z);
            }
    scanf("%d",&n);
    }
}
