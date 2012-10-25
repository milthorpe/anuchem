#!/bin/bash

PARTITION=$1
MAX_PLACES=$2
n=200000
dmax=4
p=12
ws=1

i=1
while [ $i -le $MAX_PLACES ]
do
  for j in 1 2 3
  do
# -env BG_COREDUMPONEXIT=1
    echo "mpirun -env X10RT_PGAS_COLLECTIVES=true -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/fmm $n $dmax $p $ws -verbose"
    mpirun -env X10RT_PGAS_COLLECTIVES=true -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/fmm $n $dmax $p $ws -verbose 
    echo ""
  done
  i=$(( i*2 ))
done
