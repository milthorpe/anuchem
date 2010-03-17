#ifndef ORG_GNU_GSL_H
#define ORG_GNU_GSL_H

#include "gsl/gsl_math.h"  
#include "gsl/gsl_blas.h"  
#include "gsl/gsl_eigen.h" 
#include "gsl/gsl_linalg.h"

namespace org {
   namespace gnu {   
      namespace gsl {
         class GSLWrapper {
            public:
               static int eigenSymmv(double *A, double *eval, double *evec, int aSize);
               static int solve(double *A, double *b, double *x, int aSize);
         }; 
      } // gsl
   } // gnu
} // org

#endif

