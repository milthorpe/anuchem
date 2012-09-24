N=60000
Q=10
DMAX=3
P=12
WS=1

echo "X10_NTHREADS=8 bin/fmm $N $Q $DMAX $P $WS -verbose"
X10_NTHREADS=8 bin/fmm $N $Q $DMAX $P $WS -verbose
