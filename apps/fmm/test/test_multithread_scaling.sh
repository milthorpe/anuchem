#!/bin/bash

p=6
i=1000000
ws=1
dmax=5
for  ncores in 1 2 4 8
do

    echo "X10_NTHREADS=$ncores bin/fmm   $i $dmax $p $ws -verbose"
    for  j in 1 2
    do
        X10_NTHREADS=$ncores bin/fmm   $i $dmax $p $ws -verbose
    done
done

