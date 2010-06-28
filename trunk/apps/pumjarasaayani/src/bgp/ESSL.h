#ifndef BGP_ESSL_H
#define BGP_ESSL_H

#include "essl.h"  

namespace bgp {
   namespace essl {   
         class ESSLWrapper {
            public:
               static int eigenSymmv(double *A, double *eval, double *evec, int aSize);
               static int solve(double *A, double *b, double *x, int aSize);
         }; 
   } // essl
} // bgp

#endif

