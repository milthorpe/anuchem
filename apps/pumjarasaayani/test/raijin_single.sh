#!/bin/bash
# Raijin job script: single socket (1 place per 8 core Sandy Bridge)
. ~/x10-trunk/x10.profile

test[0]=water5-cc-pVDZ
test[1]=water5-cc-pVTZ
test[2]=water5-cc-pVQZ
test[3]=water10-cc-pVDZ
test[4]=water10-cc-pVTZ
test[5]=water10-cc-pVQZ
test[6]=1D-tetrapep-6311g
test[7]=1D-tetrapep-cc-pVDZ
test[8]=1D-tetrapep-cc-pVTZ
test[9]=1D-tetrapep-cc-pVQZ
test[10]=3D-tetrapep-6311g
test[11]=3D-tetrapep-cc-pVDZ
test[12]=3D-tetrapep-cc-pVTZ
test[13]=3D-tetrapep-cc-pVQZ
test[14]=1D-octapep-6311g
test[15]=3D-octapep-6311g
test[16]=3D-octapep-cc-pVDZ
test[17]=3D-octapep-cc-pVTZ
test[18]=3D-octapep-cc-pVQZ

for i in {0..18}
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
