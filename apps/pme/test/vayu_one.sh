#!/bin/bash
# Vayu job script: test PME on one node (2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=8GB,ncpus=8
#PBS -wd
cat /proc/cpuinfo > cpuinfo.out
echo "1 CPU"
X10_NTHREADS=1 bin/pme 51396 
echo "\n1/8" 
X10_NTHREADS=8 bin/pme 51396 
X10_NTHREADS=8 bin/pme 51396 
X10_NTHREADS=8 bin/pme 51396 
echo "\n2/4" 
mpiexec -n 2 -x X10_NTHREADS=4 bin/pme 51396 
mpiexec -n 2 -x X10_NTHREADS=4 bin/pme 51396 
mpiexec -n 2 -x X10_NTHREADS=4 bin/pme 51396 
echo "\n4/2" 
mpiexec -n 4 -x X10_NTHREADS=2 bin/pme 51396 
mpiexec -n 4 -x X10_NTHREADS=2 bin/pme 51396 
mpiexec -n 4 -x X10_NTHREADS=2 bin/pme 51396 
echo "\n8/1" 
mpiexec -n 8 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396 
mpiexec -n 8 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396 
mpiexec -n 8 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pme 51396 

echo "\n16 threads" 

echo "\n1/16" 
X10_NTHREADS=16 bin/pme 51396 
X10_NTHREADS=16 bin/pme 51396 
X10_NTHREADS=16 bin/pme 51396 
echo "\n2/8" 
mpiexec -n 2 -x X10_NTHREADS=8 bin/pme 51396 
mpiexec -n 2 -x X10_NTHREADS=8 bin/pme 51396 
mpiexec -n 2 -x X10_NTHREADS=8 bin/pme 51396 
echo "\n4/4" 
mpiexec -n 4 -x X10_NTHREADS=4 bin/pme 51396 
mpiexec -n 4 -x X10_NTHREADS=4 bin/pme 51396 
mpiexec -n 4 -x X10_NTHREADS=4 bin/pme 51396 
echo "\n8/2" 
mpiexec -n 8 -x X10_NTHREADS=2 bin/pme 51396 
mpiexec -n 8 -x X10_NTHREADS=2 bin/pme 51396 
mpiexec -n 8 -x X10_NTHREADS=2 bin/pme 51396 
