module unload openmpi
module load intel-mpi
for iter in {1..10}
do
  echo "mpirun -perhost 2 -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani test/water5-cc-pVQZ.inp"
  mpirun -perhost 2 -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani test/water5-cc-pVQZ.inp
done
