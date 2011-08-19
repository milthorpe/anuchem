X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -x10rt mpi
RUNX10FLAGS := X10_NTHREADS=2 X10_STATIC_THREADS=true 
X10_TEST_SOURCE := ~/x10-trunk/x10.tests/examples/x10lib

run:    parallelLoop
	${RUNX10FLAGS} ./parallelLoop 8000

parallelLoop:	BenchmarkParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) -sourcepath ${X10_TEST_SOURCE} BenchmarkParallelLoop.x10 -o parallelLoop

clean:
	rm -rf parallelLoop *.cc *.h *.dat *.ps harness

