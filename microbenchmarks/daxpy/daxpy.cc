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
#include <algorithm>

const size_t ITERS = 1000;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

/** 
 * Simple code to perform a dot product of vectors of length N.
 * Intended as a comparison to DotProduct.x10
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: dotProd N\n");
        exit(1);
    }
    int N = atoi(argv[1]);
    double alpha = 2.5;
    double *x = new double[N]();
    double *y = new double[N]();
    for (size_t i = 0; i < N; i++) {
        x[i] = y[i] = N;
    }

    int64_t t1 = TimeInMicros();
    for (size_t iter = 0; iter < ITERS; iter++) {
        for (size_t i = 0; i < N; i++) {
            x[i] = alpha * x[i] + y[i];
        }
    }
    int64_t t2 = TimeInMicros();

    printf("C++ DAXPY for vectors length %d time %f ms\n", N, (t2-t1)/1e3/ITERS);

    delete x;
    delete y;
}
