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
i=5
while [ $i -le 30 ]
do
  run_test "X10_NTHREADS=4 bin/fmm3d 10000 60 $i 2 -verbose -compare"
  i=$(( i+5 ))
done

run_test "X10_NTHREADS=4 bin/fmm3d 51396 60 10 2 -verbose -compare"
run_test "mpiexec -n 2 -x X10_NTHREADS=2 bin/fmm3d 51396 60 10 2 -verbose -compare"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 bin/fmm3d 51396 60 10 2 -verbose -compare"
run_test "mpiexec -n 4 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/fmm3d 51396 60 10 2 -verbose -compare"
run_test "X10_NTHREADS=4 bin/fmm3d 100000 60 10 2 -verbose -compare"
run_test "X10_NTHREADS=4 bin/fmm3d 250000 64 10 2 -verbose -compare"

