X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -x10rt mpi #-WORK_STEALING=true
RUNX10FLAGS := X10_NTHREADS=4 X10_STATIC_THREADS=true 

run:    parallelLoop
	${RUNX10FLAGS} ./parallelLoop 8000

parallelLoop:	BenchmarkParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkParallelLoop.x10 -o parallelLoop

clean:
	rm -rf parallelLoop *.cc *.h *.dat *.ps harness

