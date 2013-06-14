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
#include <tbb/tbb.h>
#include <algorithm>

const size_t N = 10000;
const size_t ITERS = 1000;

int64_t TimeInMicros() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

/** 
 * Simple code to perform a parallel sum of N copies of 1.0.
 * Intended as a comparison to BenchmarkWorkerLocalReduce.x10
 */
int main(int argc, char *argv[]) {
    double sum = 0.0;
    double totalsum = 0.0;
    int64_t t1 = TimeInMicros();
    for (size_t i = 0; i < ITERS; i++) {
        sum = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0, N, 1),
            0.0,
            [](const tbb::blocked_range<size_t>& r, double value)->double
            {
                return value + (double)(r.end() - r.begin());
            },
            std::plus<double>(),
            tbb::simple_partitioner()
        );
        totalsum += sum;
    }
    int64_t t2 = TimeInMicros();

    printf("sum %f totalsum %f time %f ms\n", sum, totalsum, (t2-t1)/1e3/ITERS);
}
