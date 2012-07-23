#!/bin/bash

i=5000
dmax=3
p=12
while [ $i -le 100000 ]
do
  echo "mpiexec -n 8 ./directcpp $i"
  for j in 1 2 3
  do
      mpiexec -n 8 ./directcpp $i
  done
  i=$(( i+5000 ))
done

