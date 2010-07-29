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

#ifndef BGP_ESSL_H
#define BGP_ESSL_H

#include "essl.h"  

/**
 * ESSL wrapper header
 *
 * @author V. Ganesh
 */
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

