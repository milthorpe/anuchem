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

const size_t ITERS = 1000;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

struct Daxpy {
    double* x;
    const double *y;
    double alpha;

    Daxpy(const double alpha, double* x, const double* y) {
        this->alpha = alpha;
        this->x = x;
        this->y = y;
    }

    void operator()(const tbb::blocked_range<size_t>& range) const {
        for (int i = range.begin(); i < range.end(); ++i) {
            x[i] = alpha * x[i] + y[i];
        }
    }
};

/** 
 * C++/TBB code to perform a DAXPY of vectors of length N.
 * Intended as a comparison to Daxpy.x10
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "usage: daxpy_c N\n");
        exit(1);
    }
    printf("C++ DAXPY with TBB\n");
    int N = atoi(argv[1]);
    int num_threads =  tbb::task_scheduler_init::default_num_threads();
    if (argc > 2) {
        num_threads = atoi(argv[2]);
    }
    double alpha = 2.5;
    double *x = new double[N]();
    double *y = new double[N]();
    for (size_t i = 0; i < N; i++) {
        x[i] = y[i] = N;
    }

    tbb::task_scheduler_init sched(num_threads);

    Daxpy daxpy(alpha, x, y);

    int64_t t1 = TimeInMicros();
    for (size_t iter = 0; iter < ITERS; iter++) {
        tbb:parallel_for(tbb::blocked_range<size_t>(0, N), daxpy);
    }
    int64_t t2 = TimeInMicros();

    printf("vector size: %d num threads: %d time per iteration: %g ms\n",
        N, num_threads, (t2-t1)/1e3/ITERS);

    delete x;
    delete y;
}
