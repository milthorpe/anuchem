#!/bin/bash
# Vayu job script: test PME scaling to 256 places (8 places per 2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
module load openmpi
cat /proc/cpuinfo > cpuinfo.out
i=1
while [ $i -le 256 ]
do
  echo "\n$i places"
  for j in 1 2
  do
    echo "mpiexec -n $i -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396"
    mpiexec -n $i -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396
    echo ""
  done
  i=$(( i*2 ))
done

