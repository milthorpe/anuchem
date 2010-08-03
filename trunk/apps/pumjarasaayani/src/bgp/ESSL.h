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

