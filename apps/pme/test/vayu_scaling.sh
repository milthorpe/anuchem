#!/bin/bash
# Vayu job script: test PME scaling on eight nodes (2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=64GB,ncpus=64
#PBS -wd
module load openmpi
cat /proc/cpuinfo > cpuinfo.out
echo "1 core"
mpiexec -n 1 -x X10_NTHREADS=1 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=1 bin/pme 51396
echo "\n2 cores"
mpiexec -n 1 -x X10_NTHREADS=2 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=2 bin/pme 51396
echo "\n3 cores"
mpiexec -n 1 -x X10_NTHREADS=3 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=3 bin/pme 51396
echo "\n4 cores"
mpiexec -n 1 -x X10_NTHREADS=4 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=4 bin/pme 51396
echo "\n5 cores"
mpiexec -n 1 -x X10_NTHREADS=5 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=5 bin/pme 51396
echo "\n6 cores"
mpiexec -n 1 -x X10_NTHREADS=6 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=6 bin/pme 51396
echo "\n7 cores"
mpiexec -n 1 -x X10_NTHREADS=7 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=7 bin/pme 51396
echo "\n8 cores"
mpiexec -n 1 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 1 -x X10_NTHREADS=8 bin/pme 51396
echo ""
echo "\n16 cores"
mpiexec -n 2 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 2 -x X10_NTHREADS=8 bin/pme 51396
echo "\n24 cores"
mpiexec -n 3 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 3 -x X10_NTHREADS=8 bin/pme 51396
echo "\n32 cores"
mpiexec -n 4 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 4 -x X10_NTHREADS=8 bin/pme 51396
echo "\n40 cores"
mpiexec -n 5 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 5 -x X10_NTHREADS=8 bin/pme 51396
echo "\n48 cores"
mpiexec -n 6 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 6 -x X10_NTHREADS=8 bin/pme 51396
echo "\n56 cores"
mpiexec -n 7 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 7 -x X10_NTHREADS=8 bin/pme 51396
echo "\n64 cores"
mpiexec -n 8 -x X10_NTHREADS=8 bin/pme 51396
mpiexec -n 8 -x X10_NTHREADS=8 bin/pme 51396

