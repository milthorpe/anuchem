#!/bin/bash

run_test() {
  for j in 1 2 3
  do
    echo "$1"
    echo "$(eval $1)"
    echo ""
  done
}

i=5000
while [ $i -le 100000 ]
do
  run_test "mpiexec -n 2 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/fmm3d $i 60 21"
  i=$(( i+5000 ))
done


