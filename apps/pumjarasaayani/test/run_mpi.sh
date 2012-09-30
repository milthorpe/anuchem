cat /proc/cpuinfo > cpuinfo.out

# TODO multithreading set X10_NTHREADS=8
echo "mpiexec -n $X10_NPLACES -x X10_NTHREADS=1 bin/pumjarasaayani test/benzene-321g.inp"
mpiexec -n $X10_NPLACES -x X10_NTHREADS=1 bin/pumjarasaayani test/benzene-321g.inp
