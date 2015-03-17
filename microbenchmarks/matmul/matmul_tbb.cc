/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2015.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <tbb/tbb.h>
#include <algorithm>

const size_t ITERS = 100;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

struct MatMul {
    const double *a, *b;
    double *c;
    const size_t M, K;

    MatMul(const double *a, const double *b, double *c, size_t M, size_t K) :
        a(a), b(b), c(c), M(M), K(K)
    {}

    void operator()(const tbb::blocked_range2d<size_t>& r) const {
        for (size_t i = r.rows().begin(); i < r.rows().end(); i++){
            for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
                double x = 0.0;
                for (size_t k = 0; k < K; k++) {
                    x += a[i+k*M] * b[k+j*K];
                }
                c[i+j*M] = x;
            }
        }
    }
};

/** 
 * C++/TBB code to perform multiplication of square matrices of dimension N.
 * Intended as a comparison to MatMul.x10
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: matmul_c N\n");
        exit(1);
    }
    printf("C++ matrix multiplication with TBB\n");
    size_t N = atoi(argv[1]);
    int num_threads =  tbb::task_scheduler_init::default_num_threads();
    if (argc > 2) {
        num_threads = atoi(argv[2]);
    }
    double *a = new double[N*N]();
    double *b = new double[N*N]();
    double *c = new double[N*N]();
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            a[i+j*N] = b[i+j*N] = i+j*N;
        }
    }

    tbb::task_scheduler_init sched(num_threads);

    MatMul matmul(a, b, c, N, N);

    int64_t t1 = TimeInMicros();
    for (size_t iter = 0; iter < ITERS; iter++) {
        tbb:parallel_for(tbb::blocked_range2d<size_t>(0, N, 0, N), matmul);
    }
    int64_t t2 = TimeInMicros();

    printf("matrix size: %d num threads: %d time per iteration: %g ms\n",
        N, num_threads, (t2-t1)/1e3/ITERS);

    delete a;
    delete b;
    delete c;
}
