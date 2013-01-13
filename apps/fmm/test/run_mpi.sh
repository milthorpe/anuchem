N=1000000
DMAX=5
P=6
WS=1

cat $PBS_NODEFILE > nodefile.$X10_NPLACES.out

cat /proc/cpuinfo > cpuinfo.out

echo "mpiexec -n $X10_NPLACES -bind-to-socket -cpus-per-proc $X10_NTHREADS -x X10_NTHREADS=$X10_NTHREADS bin/fmm $N $DMAX $P $WS -verbose"
mpiexec -n $X10_NPLACES -bind-to-socket -cpus-per-proc $X10_NTHREADS -x X10_NTHREADS=$X10_NTHREADS bin/fmm $N $DMAX $P $WS -verbose
