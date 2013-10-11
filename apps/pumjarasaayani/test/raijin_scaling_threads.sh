#!/bin/bash
# Raijin test script: test strong scaling with threads on a single 16-core node
#module load ipm
#module load hpctoolkit
i=1
while [ $i -le 18 ]
do
  export X10_NPLACES=1
  export X10_NTHREADS=$i
  export NCPUS=$((X10_NPLACES*16))
  MEM=$((NCPUS*2))GB
  echo $NCPUS $MEM
  qsub -l wd -lwalltime=00:15:00 -lncpus=$NCPUS -l mem=$MEM -V test/run_ten.sh
  echo ""
  i=$(( i+1 ))
done

