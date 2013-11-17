#!/bin/bash
# Raijin job script: single socket (1 place per 8 core Sandy Bridge)
. ~/x10-trunk/x10.profile

test[0]=1D-tetrapep-6311g
test[1]=1D-tetrapep-cc-pVDZ
test[2]=1D-tetrapep-cc-pVTZ
test[3]=3D-tetrapep-cc-pVDZ
test[4]=3D-tetrapep-cc-pVTZ

for i in {0..4}
do
    export TESTFILE=test/${test[i]}.inp

    export X10_NPLACES=1
    export X10_NTHREADS=1
    export NCPUS=$((X10_NPLACES*8))
    MEM=$((NCPUS*2))GB # could be 3GB/CPU, but go easy
    echo $NCPUS $MEM
    qsub -l wd -lwalltime=06:00:00 -lncpus=$NCPUS -l mem=$MEM -V test/run_single.sh
    echo ""
done
