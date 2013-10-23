. ~/x10-trunk/x10.profile
module unload openmpi
module load intel-mpi
uniq < $PBS_NODEFILE > hosts.txt

N=1000000
DMAX=5
P=6
WS=1

echo "mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/fmm $N $DMAX $P $WS -verbose"
mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/fmm $N $DMAX $P $WS -verbose
