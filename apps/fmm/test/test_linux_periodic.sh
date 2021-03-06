#!/bin/bash

run_test() {
  for j in 1 2 3
  do
    echo "$1"
    echo "$(eval $1)"
    echo ""
  done
}

unset X10_STATIC_THREADS
run_test "X10_NTHREADS=4 bin/periodicFmm3d 51396 60 10 10 -verbose"
run_test "mpiexec -n 2 -x X10_NTHREADS=2 bin/periodicFmm3d 51396 60 10 10 -verbose"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 bin/periodicFmm3d 51396 60 10 10 -verbose"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/periodicFmm3d 51396 60 10 10 -verbose"




