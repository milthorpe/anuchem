#!/bin/bash

p=20
start3=35000
i=5000
dmax=2
ws=1
ncores=8
while [ $i -le 100000 ]
do
  echo "X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws -verbose -compare -forces"
  X10_NTHREADS=1 X10_STATIC_THREADS=true bin/fmm $i 150 $dmax $p $ws -verbose -compare -forces
  echo "mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws"
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws #-verbose -compare
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws
  mpiexec -n $ncores -x X10_NTHREADS=1 bin/fmm $i 150 $dmax $p $ws
  i=$(( i+5000 ))
  if [ $i -ge $start3 ]; then dmax=3; fi
done

