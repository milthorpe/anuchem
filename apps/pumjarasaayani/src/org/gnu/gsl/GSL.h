/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */

#ifndef ORG_GNU_GSL_H
#define ORG_GNU_GSL_H

#include <x10/array/Array.h>

#include "gsl/gsl_math.h"  
#include "gsl/gsl_blas.h"  
#include "gsl/gsl_eigen.h" 
#include "gsl/gsl_linalg.h"

/**
 * Header for GSL wrapper
 *
 * @author V. Ganesh
 */
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

