#!/bin/bash
# Vayu job script: test FMM scaling to 256 places (1 place per 4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
#module load ipm
#module load hpctoolkit
i=1
while [ $i -le 256 ]
do
  export X10_NPLACES=$i
  export X10_NTHREADS=4
  export NCPUS=$((X10_NPLACES*X10_NTHREADS))
  VMEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
  echo $NCPUS $VMEM
  qsub -wd -lwalltime=00:02:00 -lncpus=$NCPUS -lvmem=$VMEM -V test/run_mpi.sh
  echo ""
  i=$(( i*2 ))
done

