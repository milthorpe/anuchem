. ~/x10-trunk/x10.profile
module unload openmpi
module load intel-mpi
export HOSTS=${PBS_JOBID}_hosts.txt
uniq < $PBS_NODEFILE > $HOSTS
for iter in {1..10}
do
  echo "mpirun -binding domain=socket -perhost 2 -hostfile $HOSTS -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani test/water10-cc-pVQZ.inp"
  mpirun -binding domain=socket -perhost 2 -hostfile $HOSTS -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani test/water10-cc-pVQZ.inp
done
