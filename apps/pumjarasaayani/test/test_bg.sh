#!/bin/bash

PARTITION=$1
MAX_PLACES=$2

i=1
while [ $i -le 128 ]
do
  for j in 1 
  do
    echo "mpirun -env X10RT_PGAS_COLLECTIVES=true -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/pumjarasaayani test/benzene-321g.inp"
    mpirun -env X10RT_PGAS_COLLECTIVES=true -env X10_NTHREADS=1 -env X10_STATIC_THREADS=true -noallocate -nofree -partition $PARTITION -mode VN -np $i bin/pumjarasaayani test/benzene-321g.inp 
    echo ""
  done
  i=$(( i*2 ))
done
