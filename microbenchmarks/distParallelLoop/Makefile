X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -x10rt mpi
RUNX10FLAGS := X10_NTHREADS=2 X10_STATIC_THREADS=true 
X10_TEST_SOURCE := ~/x10-trunk/x10.tests/examples/x10lib

run:    distParallelLoop
	${RUNX10FLAGS} ./distParallelLoop 8000

distParallelLoop:	BenchmarkDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) -sourcepath ${X10_TEST_SOURCE} BenchmarkDistParallelLoop.x10 -o distParallelLoop

plot:	localityDistParallelLoop
	${RUNX10FLAGS} ./localityDistParallelLoop 8000 -print
	gnuplot < scatter.plot

localityParallelLoop:	LocalityDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) -sourcepath ${X10_TEST_SOURCE} LocalityDistParallelLoop.x10 -o localityDistParallelLoop

clean:
	rm -rf distParallelLoop localityParallelLoop *.cc *.h *.dat *.ps harness

