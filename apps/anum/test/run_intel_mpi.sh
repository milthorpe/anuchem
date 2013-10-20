. ~/x10-trunk/x10.profile
module unload openmpi
module load intel-mpi
uniq < $PBS_NODEFILE > hosts.txt

echo "mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/anumm test/two_aminos.inp"
mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/anumm test/two_aminos.inp
