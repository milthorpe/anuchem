#!/bin/bash

p=15
i=10000
dmax=3
ncores=8
ws=1
while [ $i -le 100000 ]
do
  echo "X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws"
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws
  X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws
  i=$(( i+10000 ))
done

