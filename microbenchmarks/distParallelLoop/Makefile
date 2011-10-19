X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -WORK_STEALING=true -x10rt mpi
RUNX10FLAGS := X10_NTHREADS=4 X10_STATIC_THREADS=true 

run:    distParallelLoop
	${RUNX10FLAGS} ./distParallelLoop 8000

distParallelLoop:	BenchmarkDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkDistParallelLoop.x10 -o distParallelLoop

plot:	localityDistParallelLoop
	${RUNX10FLAGS} ./localityDistParallelLoop 8000 -print
	gnuplot < scatter.plot

localityDistParallelLoop:	LocalityDistParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) LocalityDistParallelLoop.x10 -o localityDistParallelLoop

clean:
	rm -rf distParallelLoop localityDistParallelLoop *.cc *.h *.dat *.ps harness

