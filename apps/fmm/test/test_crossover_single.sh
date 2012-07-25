#!/bin/bash

i=20000
dmax=3
p=15
while [ $i -le 100000 ]
do
  #echo "X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2 -verbose -compare -forces"
  #X10_NTHREADS=1 X10_STATIC_THREADS=true bin/fmm $i 150 $dmax $p 2 -verbose -compare -forces
  echo "mpiexec -n 8 -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2"
  mpiexec -n 8 -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2 # -verbose -compare
  mpiexec -n 8 -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  mpiexec -n 8 -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  i=$(( i+5000 ))
done
