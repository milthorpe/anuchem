N=60000
DMAX=3
P=12
WS=1

echo "X10_NTHREADS=8 bin/fmm $N $DMAX $P $WS -verbose"
X10_NTHREADS=8 bin/fmm $N $DMAX $P $WS -verbose
