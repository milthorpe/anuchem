#!/bin/bash

PARTITION=$1
MAX_PLACES=$2
X10_NTHREADS=16

N=1000000
DMAX=5
P=6
WS=1

i=1
while [ $i -le $MAX_PLACES ]
do
  echo "runjob --block $PARTITION --ranks-per-node 1 --args $N $DMAX $P $WS --args '-verbose' --exe bin/fmm--envs X10_NTHREADS=$X10_NTHREADS --envs X10_STATIC_THREADS=false --np=$i"
  runjob --block $PARTITION --ranks-per-node 1 --args $N $DMAX $P $WS --args '-verbose' --exe bin/fmm --envs X10_NTHREADS=$X10_NTHREADS --envs X10_STATIC_THREADS=false --np=$i
  echo ""
  i=$(( i*2 ))
done
