#!/bin/bash

p=15
start3=35000
i=5000
dmax=2
ncores=8
while [ $i -le 100000 ]
do
  #echo "X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2 -verbose -compare -forces"
  #X10_NTHREADS=1 X10_STATIC_THREADS=true bin/fmm $i 150 $dmax $p 2 -verbose -compare -forces
  echo "mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2"
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2 -verbose -compare
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p 2
  i=$(( i+5000 ))
  if [ $i -ge $start3 ]; then dmax=3; fi
done

