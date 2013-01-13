cat $PBS_NODEFILE > nodefile.$X10_NPLACES.out

cat /proc/cpuinfo > cpuinfo.out

echo "mpiexec -n $X10_NPLACES -bind-to-socket -cpus-per-proc $X10_NTHREADS -x X10_NTHREADS=$X10_NTHREADS bin/anumm test/two_aminos.inp"
mpiexec -n $X10_NPLACES -bind-to-socket -cpus-per-proc $X10_NTHREADS -x X10_NTHREADS=$X10_NTHREADS bin/anumm test/two_aminos.inp
