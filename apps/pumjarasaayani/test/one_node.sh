#!/bin/bash
# Vayu job script: test Pumja Rasayaani on one node (2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=8GB,ncpus=8
#PBS -wd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/gsl/lib
cat /proc/cpuinfo > cpuinfo.out
echo "1 CPU" > test.out
X10_NTHREADS=1 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n1/8" >> test.out
X10_NTHREADS=8 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n2/4" >> test.out
mpiexec -n 2 -x X10_NTHREADS=4 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n4/2" >> test.out
mpiexec -n 4 -x X10_NTHREADS=2 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n8/1" >> test.out
mpiexec -n 8 -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pumjarasaayani test/benzene-321.inp 5 >> test.out

echo "\n16 threads" >> test.out

echo "\n1/16" >> test.out
X10_NTHREADS=16 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n2/8" >> test.out
mpiexec -n 2 -x X10_NTHREADS=8 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n4/4" >> test.out
mpiexec -n 4 -x X10_NTHREADS=4 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
echo "\n8/2" >> test.out
mpiexec -n 8 -x X10_NTHREADS=2 bin/pumjarasaayani test/benzene-321.inp 5 >> test.out
