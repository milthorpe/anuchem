/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

#ifndef ORG_GNU_GSL_H
#define ORG_GNU_GSL_H

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

