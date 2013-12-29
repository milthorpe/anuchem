. ~/x10-trunk/x10.profile
module unload openmpi
module load intel-mpi
export HOSTS=${PBS_JOBID}_hosts.txt
uniq < $PBS_NODEFILE > $HOSTS
export TESTFILE=test/3D-tetrapep-cc-pVQZ.inp
for iter in {1..10}
do
  echo "mpirun -binding domain=socket -perhost 2 -hostfile $HOSTS -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani $TESTFILE"
  mpirun -binding domain=socket -perhost 2 -hostfile $HOSTS -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pumjarasaayani $TESTFILE
done
