#!/bin/bash
# Raijin job script: test PR scaling to 8 places (1 place per 8 core Sandy Bridge)
. ~/x10-trunk/x10.profile
i=1
while [ $i -le 128 ]
do
  export X10_NPLACES=$i
  export X10_NTHREADS=8
  export NCORES=$((X10_NPLACES*8))
  MEM=$((NCORES*2))GB 
  echo $NCORES $MEM
  qsub -l wd -lwalltime=00:45:00 -lncpus=$NCORES -l mem=$MEM -V test/run_intel_mpi.sh
  echo ""
  i=$(( i*2 ))
done

