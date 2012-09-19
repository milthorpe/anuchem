#!/bin/bash

p=12
i=10000
dmax=3
ncores=8
while [ $i -le 100000 ]
do
  echo "X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2"
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  i=$(( i+10000 ))
done

