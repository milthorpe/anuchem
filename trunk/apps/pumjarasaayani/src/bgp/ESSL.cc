#include "ESSL.h"

namespace bgp {
   namespace essl {   
         int ESSLWrapper::eigenSymmv(double *A, double *evec, double *eval, int aSize) {
             int N = aSize*aSize;
             double dummy;
             double *Acpy = (double *) malloc(sizeof(double)*N);
             double *aux  = (double *) malloc(sizeof(double)*aSize*2); 

             // make copy of A the original matrix is destroyed
             for(int i=0; i<N; i++) Acpy[i] = A[i];

             // compute eigen values, vectors
             dgeev(1, Acpy, aSize, eval, evec, aSize, &dummy, aSize, aux,  aSize*2); 

             // TODO: check correctness
             // sort it
             int *indx = (int *) malloc(sizeof(int)*aSize);
             dsortx(eval, 2, aSize, indx);
             for(int i=0; i<aSize; i++) {
                 int orgIdx = indx[i];

                 // swap i and orgIdx in evec
                 for(int j=0; j<aSize; j++) {
                     dummy = evec[j][i];
                     evec[j][i] = evec[j][orgIdx]; 
                     evec[j][orgIdx] = dummy;
                 } // end for
             } // end for
             free(indx);

             free(Acpy);
             free(aux);

             return 0;
         }

         int ESSLWrapper::solve(double *A, double *b, double *x, int aSize) {
             // TODO:
             int *ipvt = (int *) malloc(sizeof(int)*aSize);
             double *Acpy = (double *) malloc(sizeof(double)*N);
             int info;

             for(int i=0; i<N; i++) Acpy[i] = A[i];

             dgetrf(aSize, aSize, A, aSize, ipvt, &info);
             dgetrs('N', aSize, aSize, aSize, aSize, ipvt, BX, aSize, &info);

             free(ipvt);
             free(Acpy);

             return info;
         }
   } // essl
} // bgp


