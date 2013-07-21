#!/bin/bash
# Vayu job script: test PME scaling to 256 places (8 places per 2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
#module load ipm
i=1
while [ $i -le 256 ]
do
  export X10_NPLACES=$i
  export VMEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
  qsub -wd -lwalltime=00:10:00 -lncpus=$i -lvmem=$VMEM -V test/run_three.sh
  echo ""
  i=$(( i*2 ))
done

