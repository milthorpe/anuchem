. ~/x10-trunk/x10.profile
module unload openmpi
module load intel-mpi
uniq < $PBS_NODEFILE > hosts.txt

MOLFILE=test/gromacs/box100k.gro
BETA=0.35
CUTOFF=10.0
GRIDSIZE=96

echo "mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pme $MOLFILE $BETA $CUTOFF $GRIDSIZE"
mpirun -binding domain=socket -perhost 2 -hostfile hosts.txt -n $X10_NPLACES -env X10_NTHREADS $X10_NTHREADS -env OMP_NUM_THREADS $X10_NTHREADS bin/pme $MOLFILE $BETA $CUTOFF $GRIDSIZE
