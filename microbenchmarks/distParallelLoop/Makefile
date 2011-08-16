X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -x10rt mpi
RUNX10FLAGS := X10_NTHREADS=4 X10_STATIC_THREADS=true 
X10_TEST_SOURCE := ~/x10-trunk/x10.tests/examples/x10lib

plot:	parallelLoop
	${RUNX10FLAGS} ./parallelLoop
	gnuplot < scatter.plot

parallelLoop:	BenchmarkDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) -sourcepath ${X10_TEST_SOURCE} BenchmarkDistParallelLoop.x10 -o parallelLoop

clean:
	rm -rf parallelLoop *.cc *.h *.dat *.ps harness

