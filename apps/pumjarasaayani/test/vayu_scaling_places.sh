#!/bin/bash
# Vayu job script: test PR scaling to 8 places (1 place per 4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
. ~/x10-trunk/x10.profile
i=1
while [ $i -le 8 ]
do
  export X10_NPLACES=$i
  export X10_NTHREADS=4
  export NCPUS=$((X10_NPLACES*4))
  VMEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
  echo $NCPUS $VMEM
  qsub -wd -lwalltime=00:15:00 -lncpus=$NCPUS -lvmem=$VMEM -V test/run_intel_mpi.sh
  echo ""
  i=$(( i*2 ))
done

