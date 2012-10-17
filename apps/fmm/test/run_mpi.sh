N=200000
Q=10
DMAX=4
P=12
WS=1

cat /proc/cpuinfo > cpuinfo.out

echo "mpiexec -n $X10_NPLACES -x X10_NTHREADS=8 bin/fmm $N $Q $DMAX $P $WS -verbose"
mpiexec -n $X10_NPLACES -x X10_NTHREADS=8 bin/fmm $N $Q $DMAX $P $WS -verbose
