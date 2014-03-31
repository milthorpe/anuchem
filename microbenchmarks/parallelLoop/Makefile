X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -VERBOSE_CHECKS #-WORK_STEALING=true
RUNX10ENV := X10_NTHREADS=4

run:    parallelLoop
	${RUNX10ENV} ./parallelLoop 8000

parallelLoop:	BenchmarkParallelLoop.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkParallelLoop.x10 -o parallelLoop

.PHONY: clean
clean:
	rm -rf parallelLoop *.cc *.h

