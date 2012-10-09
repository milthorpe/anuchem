#!/bin/bash
# Vayu job script: test FMM scaling to 64 places (1 place per 2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
#module load ipm
module load hpctoolkit
i=32
while [ $i -le 32 ]
do
  export X10_NPLACES=$i
  export NCPUS=$((i*8))
  VMEM=$((NCPUS*3))GB
  echo $VMEM
  qsub -wd -lwalltime=00:05:00 -lncpus=$NCPUS -lvmem=$VMEM -V test/run_mpi.sh
  echo ""
  i=$(( i*2 ))
done

