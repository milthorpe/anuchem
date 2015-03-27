X10FLAG := -O -NO_CHECKS

# enable CPU profiling with google-perftools
PROFILE ?=
ifdef PROFILE
  X10FLAG += -gpt
endif

NPROC ?= $(shell nproc)

CXX 	 = g++
CXXFLAGS = -g -O3 -fopenmp
LDFLAGS = -lpthread

SIZE=500

run:    run_x10_single run_omp_single run_tbb_single run_x10_multi run_omp_multi run_tbb_multi

run_x10_single: matmul
	X10_NTHREADS=1 ./matmul $(SIZE)

run_x10_multi: matmul
	X10_NTHREADS=$(NPROC) ./matmul $(SIZE)

run_omp_single: matmul_omp
	OMP_NUM_THREADS=1 ./matmul_omp $(SIZE)

run_omp_multi: matmul_omp
	OMP_NUM_THREADS=$(NPROC) ./matmul_omp $(SIZE)

run_tbb_single: matmul_tbb
	./matmul_tbb $(SIZE) 1

run_tbb_multi: matmul_tbb
	./matmul_tbb $(SIZE)

matmul:	MatMul.x10
	x10c++ $(X10FLAG) -post '# $(CXXFLAGS)' -report postcompile=1 MatMul.x10 -o matmul

matmul_omp:	matmul_omp.cc
		$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

matmul_tbb:	matmul_tbb.cc
		$(CXX) $(CXXFLAGS) -I$(TBB_INCLUDE) -o $@ $< $(LDFLAGS) -L$(TBB_LIB) -ltbb


clean:
	rm -rf matmul matmul_omp matmul_tbb

