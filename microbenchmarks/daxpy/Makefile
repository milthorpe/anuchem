X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS

# enable CPU profiling with google-perftools
PROFILE ?=
ifdef PROFILE
  X10FLAG += -gpt
endif

NPROC ?= $(shell nproc)

CXX 	 = g++
CXXFLAGS = -O3 -fopenmp
LDFLAGS = -lpthread

SIZE=1000000

run:    run_x10_single run_omp_single run_tbb_single run_x10_multi run_omp_multi run_tbb_multi

run_x10_single: daxpy
	X10_NTHREADS=1 ./daxpy $(SIZE)

run_x10_multi: daxpy
	X10_NTHREADS=$(NPROC) ./daxpy $(SIZE)

run_omp_single: daxpy_omp
	OMP_NUM_THREADS=1 ./daxpy_omp $(SIZE)

run_omp_multi: daxpy_omp
	OMP_NUM_THREADS=$(NPROC) ./daxpy_omp $(SIZE)

run_tbb_single: daxpy_tbb
	./daxpy_tbb $(SIZE) 1

run_tbb_multi: daxpy_tbb
	./daxpy_tbb $(SIZE)

daxpy:	Daxpy.x10
	$(X10CPP) $(X10CPPFLAGS) -post '# $(CXXFLAGS)' -report postcompile=1 $< -o daxpy

daxpy_omp:	daxpy_omp.cc
		$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

daxpy_tbb:	daxpy_tbb.cc
		$(CXX) $(CXXFLAGS) -I$(TBB_INCLUDE) -o $@ $< $(LDFLAGS) -L$(TBB_LIB) -ltbb

clean:
	rm -rf daxpy daxpy_omp daxpy_tbb

