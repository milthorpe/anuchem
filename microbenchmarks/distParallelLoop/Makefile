X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -VERBOSE_CHECKS -x10rt mpi
RUNX10ENV := X10_NTHREADS=4

run:    distParallelLoop
	${RUNX10ENV} ./distParallelLoop 8000

distParallelLoop:	BenchmarkDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkDistParallelLoop.x10 -o distParallelLoop

plot:	localityDistParallelLoop
	${RUNX10ENV} ./localityDistParallelLoop 8000 -print
	gnuplot < scatter.plot

localityDistParallelLoop:	LocalityDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) LocalityDistParallelLoop.x10 -o localityDistParallelLoop

.PHONY: clean
clean:
	rm -rf distParallelLoop localityDistParallelLoop *.cc *.h *.dat *.ps

