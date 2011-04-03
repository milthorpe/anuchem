#!/bin/bash

run_test() {
  for j in 1 2 3
  do
    echo "$1"
    echo "$(eval $1)"
    echo ""
  done
}

export X10_NTHREADS=4
unset X10_STATIC_THREADS

run_test "bin/pme 30000 0.35 10.0 40"
run_test "bin/pme 30000 0.35 10.0 64"
run_test "bin/pme 30000 0.35 10.0 80"
run_test "bin/pme 51396 0.35 10.0 64"
run_test "mpiexec -n 2 -x X10_NTHREADS=2 bin/pme 51396 0.35 10.0 64"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 bin/pme 51396 0.35 10.0 64"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396 0.35 10.0 64"
run_test "bin/pme 100000 0.35 10.0 64"

