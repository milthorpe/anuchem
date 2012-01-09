#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "bessel4.c"

//#define PI 3.14159265358979323846
#define sqr(x) ((x)*(x))
#define SQRT2 1.4142135623730951

typedef struct {int x,y,z;} Point;
int delta[3][3]={{1,0,0},{0,1,0},{0,0,1}};

namespace au {
    namespace edu {
        namespace anu {
            namespace qm {
                namespace ro {
                    class Integral_Pack {
                        public:
                            Integral_Pack(int N);
                            int Genclass(int a, int b, double *A, double *B, double *zetaA, double *zetaB, double *conA, double *conB, int dconA, int dconB, int L);

                        private:
                            int N;
                            // BRA
                            #define MAX_BRA_L 4
                            #define MAX_TOTAL_BRA_L (MAX_BRA_L+1)*(MAX_BRA_L+2)*(MAX_BRA_L+3)/6
                            int map3[MAX_BRA_L+1][MAX_BRA_L+1][MAX_BRA_L+1];
                            Point inverseMap3[MAX_TOTAL_BRA_L];
                            int buildMap[MAX_TOTAL_BRA_L];
                            int totalBraL[MAX_BRA_L+1],noOfBra[MAX_BRA_L+1];
                            Point *HRRMAP[MAX_BRA_L+1][MAX_BRA_L+1];

                            // KET
                            #define MAX_KET_L 100
                            #define MAX_KET_LM (MAX_KET_L+1)*(MAX_KET_L+1)
                            #define lm2k(l,m) ((l)*(l)+l+m)
                            double cxminus[MAX_KET_LM],cxplus[MAX_KET_LM],cyminus[MAX_KET_LM],cyplus[MAX_KET_LM],cz[MAX_KET_LM];

                            double *lambda, *q;

                            void initialize();
                            int initializeCoulomb(int N);
                            int GenJ(double *B, double x, int L);
                            int GenY(double *Y, double X, double phi, int L);

                    };
                }
            }
        }
    }
}

