/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014-2015.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>

const size_t ITERS = 100;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

/** 
 * C++/OpenMP code to perform multiplication of square matrices of dimension N.
 * Intended as a comparison to MatMul.x10
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: matmul_c N\n");
        exit(1);
    }
    printf("C++ matrix multiplication with OpenMP\n");
    int N = atoi(argv[1]);

    double *a = new double[N*N]();
    double *b = new double[N*N]();
    double *c = new double[N*N]();

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            a[i+j*N] = b[i+j*N] = i+j*N;
        }
    }

    int64_t t1 = TimeInMicros();
    for (size_t iter = 0; iter < ITERS; iter++) {
        #pragma omp for schedule(static) collapse(2)
        for (size_t j = 0; j < N; j++) {
            for (size_t i = 0; i < N; i++) {
                double x = 0.0;
                for (size_t k = 0; k < N; k++) {
                    x += a[i+k*N] * b[k+j*N];
                }
                c[i+j*N] = x;
            }
        }
    }
    int64_t t2 = TimeInMicros();

    printf("matrix size: %d num threads: %d time per iteration: %g ms\n",
        N, omp_get_max_threads(), (t2-t1)/1e3/ITERS);

    delete a;
    delete b;
    delete c;
}
