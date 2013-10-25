#!/bin/bash
# Raijin job script: test PR scaling to 8 places (1 place per 8 core Sandy Bridge)
. ~/x10-trunk/x10.profile
i=64
while [ $i -le 128 ]
do
  export X10_NPLACES=$i
  export X10_NTHREADS=8
  export NCPUS=$((X10_NPLACES*8))
  MEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
  echo $NCPUS $MEM
  qsub -l wd -lwalltime=00:10:00 -lncpus=$NCPUS -l mem=$MEM -V test/run_intel_mpi.sh
  echo ""
  i=$(( i*2 ))
done

