X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS

CXX 	 = g++
CXXFLAGS = -g -O2 -fopenmp -ftree-vectorize #-ftree-vectorizer-verbose=6 
LDFLAGS = 

run:    matmul matmul_c
	X10_NTHREADS=4 ./matmul 500; OMP_NUM_THREADS=4 ./matmul_c 500

matmul:	MatMul.x10
	$(X10CPP) $(X10CPPFLAGS) -post '# $(CXXFLAGS)' -report postcompile=1 MatMul.x10 -o matmul

matmul_c:	matmul.cc
		$(CXX) $(CXXFLAGS) -o $@ matmul.cc $(LGFLAGS)

clean:
	rm -rf matmul matmul_c

