X10FLAG := -O -NO_CHECKS

# enable CPU profiling with google-perftools
PROFILE ?=
ifdef PROFILE
  X10FLAG += -gpt
endif

NPROC ?= $(shell nproc)

CXX 	 = g++
CXXFLAGS = -g -O2 -fopenmp -ftree-vectorize #-ftree-vectorizer-verbose=6 
LDFLAGS = 

run:    run_x10_single run_cpp_single run_x10_multi run_cpp_multi

run_x10_single: matmul
	X10_NTHREADS=1 ./matmul 500

run_x10_multi: matmul
	X10_NTHREADS=$(NPROC) ./matmul 500

run_cpp_single: matmul_c
	OMP_NUM_THREADS=1 ./matmul_c 500

run_cpp_multi: matmul_c
	OMP_NUM_THREADS=$(NPROC) ./matmul_c 500

matmul:	MatMul.x10
	x10c++ $(X10FLAG) -post '# $(CXXFLAGS)' -report postcompile=1 MatMul.x10 -o matmul

matmul_c:	matmul.cc
		$(CXX) $(CXXFLAGS) -o $@ matmul.cc $(LGFLAGS)

clean:
	rm -rf matmul matmul_c

