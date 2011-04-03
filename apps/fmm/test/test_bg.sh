#!/bin/bash

run_test() {
  for j in 1 2 3
  do
    echo "$1"
    echo "$(eval $1)"
    echo ""
  done
}

PARTITION=$1
MAX_PLACES=$2

i=1
while [ $i -le 256 ]
do
  run_test "mpirun -env X10RT_PGAS_COLLECTIVES=true -env BG_COREDUMPONEXIT=1 -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/periodicFmm3d 51396 120"
  i=$(( i*2 ))
done
