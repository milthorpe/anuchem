/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2013.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>

const size_t ITERS = 1000;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

/** 
 * C++/OpenMP code to perform a DAXPY of vectors of length N.
 * Intended as a comparison to Daxpy.x10
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: daxpy_omp N\n");
        exit(1);
    }
    printf("C++ DAXPY with OpenMP\n");
    int N = atoi(argv[1]);
    double alpha = 2.5;
    double *x = new double[N]();
    double *y = new double[N]();
    for (size_t i = 0; i < N; i++) {
        x[i] = y[i] = N;
    }

    int64_t t1 = TimeInMicros();
    for (size_t iter = 0; iter < ITERS; iter++) {
        #pragma omp parallel for private(i) schedule(static)
        for (size_t i = 0; i < N; i++) {
            x[i] = alpha * x[i] + y[i];
        }
    }
    int64_t t2 = TimeInMicros();

    printf("vector size: %d OMP_NUM_THREADS: %d time per iteration: %g ms\n",
        N, omp_get_max_threads(), (t2-t1)/1e3/ITERS);

    delete x;
    delete y;
}
