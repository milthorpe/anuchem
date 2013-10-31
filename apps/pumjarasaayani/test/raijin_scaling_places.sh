#!/bin/bash
# Raijin job script: test PR scaling to 8 places (1 place per 8 core Sandy Bridge)
. ~/x10-trunk/x10.profile
i=1
while [ $i -le 16 ]
do
  export X10_NPLACES=$i
  export X10_NTHREADS=8
  export NCPUS=$((X10_NPLACES*8))
  MEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
  echo $NCPUS $MEM
  qsub -l wd -lwalltime=00:15:00 -lncpus=$NCPUS -l mem=$MEM -V test/run_intel_mpi.sh
  echo ""
  i=$(( i*2 ))
done

