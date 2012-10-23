N=200000
DMAX=4
P=12
WS=1

cat /proc/cpuinfo > cpuinfo.out

echo "mpiexec -n $X10_NPLACES -x X10_NTHREADS=8 bin/fmm $N $DMAX $P $WS -verbose"
mpiexec -n $X10_NPLACES -x X10_NTHREADS=8 bin/fmm $N $DMAX $P $WS -verbose
