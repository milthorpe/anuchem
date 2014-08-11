X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS

NPROC := $(shell nproc)

CXX 	 = g++
CXXFLAGS = -O2 -fopenmp -ftree-vectorize #-ftree-vectorizer-verbose=6 
LDFLAGS = 

run:    run_x10_single run_cpp_single run_x10_multi run_cpp_multi

run_x10_single: daxpy
	X10_NTHREADS=1 ./daxpy 1000000

run_x10_multi: daxpy
	X10_NTHREADS=$(NPROC) ./daxpy 1000000

run_cpp_single: daxpy_c
	OMP_NUM_THREADS=1 ./daxpy_c 1000000

run_cpp_multi: daxpy_c
	OMP_NUM_THREADS=$(NPROC) ./daxpy_c 1000000

daxpy:	Daxpy.x10
	$(X10CPP) $(X10CPPFLAGS) -post '# $(CXXFLAGS)' -report postcompile=1 Daxpy.x10 -o daxpy

daxpy_c:	daxpy.cc
		$(CXX) $(CXXFLAGS) -o $@ daxpy.cc $(LGFLAGS)

clean:
	rm -rf daxpy daxpy_c

