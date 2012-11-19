#!/bin/bash

i=100000
ws=1
ncores=8
dmax=3
for  p in 2 3 4 5 6 7 8 9 10
do
    echo "X10_NTHREADS=$ncores bin/fmm   $i $dmax $p $ws"
    for  j in 1 2 3
    do
        X10_NTHREADS=$ncores bin/fmm   $i $dmax $p $ws
    done
    echo ""
done
#while [ $i -le 1000000 ]
#do
#  X10_NTHREADS=$ncores bin/fmm $i 150 $dmax $p $ws
#  i=$(( i+100000 ))
#done

