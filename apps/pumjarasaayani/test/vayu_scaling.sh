#!/bin/bash
# Vayu job script: test Pumja Rasaayani scaling to N places (1 place per 2x4 core Nehalem)
#PBS -P y42 
#PBS -q normal 
#PBS -l walltime=01:00:00,vmem=256GB,ncpus=256
#PBS -wd
#uncomment following line for IPM (MPI profiling)
#module load ipm
#uncomment following line for HPCToolkit (multi-process CPU profiling and tracing)
#module load hpctoolkit
i=1
while [ $i -le 1 ]
do
  export X10_NPLACES=$i
  export NCPUS=1 #$((i*8)) TODO multithreaded
  VMEM=$((NCPUS*3))GB
  echo $VMEM
  qsub -wd -lwalltime=00:05:00 -lncpus=$NCPUS -lvmem=$VMEM -V test/run_mpi.sh
  echo ""
  i=$(( i*2 ))
done

