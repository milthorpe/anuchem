#!/bin/bash
# Vayu job script: test Pumja Rasayaani on one node (2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=64GB,ncpus=64
#PBS -wd
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/gsl/lib
cat /proc/cpuinfo > cpuinfo.out

i=1
while [ $i -le 64 ]
do
    echo "mpirun -env X10RT_PGAS_COLLECTIVES=true -env BG_COREDUMPONEXIT=1 -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/pumjarasaayani test/benzene-321.inp 5"
    mpiexec -n $i -x X10_NTHREADS=1 -x X10_STATIC_THREADS=true bin/pumjarasaayani test/benzene-321.inp 5 
    echo ""
  i=$(( i*2 ))
done

