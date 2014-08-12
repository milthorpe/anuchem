X10FLAG := -O -NO_CHECKS #-WORK_STEALING=true
RUNX10FLAGS := X10_NTHREADS=4

# enable CPU profiling with google-perftools
PROFILE ?=
ifdef PROFILE
  X10FLAG += -gpt
endif

CXX 	 = g++
TBB_CXXFLAGS = -std=c++0x -O2 -I /home/josh/tbb/include/
TBB_LDFLAGS = -L /home/josh/tbb/build/linux_intel64_gcc_cc4.6_libc2.15_kernel3.2.0_debug -ltbb

run:    workerLocalReduce
	${RUNX10FLAGS} ./workerLocalReduce 10000

workerLocalReduce:	BenchmarkWorkerLocalReduce.x10
	x10c++ $(X10FLAG) BenchmarkWorkerLocalReduce.x10 -o workerLocalReduce

parallel_sum:	parallel_sum.cc
		$(CXX) $(TBB_CXXFLAGS) -o $@ parallel_sum.cc $(TBB_LDFLAGS)

clean:
	rm -rf parallel_sum workerLocalReduce BenchmarkWorkerLocalReduce.cc BenchmarkWorkerLocalReduce.h *.dat *.ps harness

