X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS

CXX 	 = g++
CXXFLAGS = -O2 -ftree-vectorizer-verbose=6 -ftree-vectorize 
LDFLAGS = 

run:    daxpy daxpy_c
	X10_NTHREADS=1 ./daxpy 1000000; ./daxpy_c 1000000

daxpy:	Daxpy.x10
	$(X10CPP) $(X10CPPFLAGS) -post '# $(CXXFLAGS)' -report postcompile=1 Daxpy.x10 -o daxpy

daxpy_c:	daxpy.cc
		$(CXX) $(CXXFLAGS) -o $@ daxpy.cc $(LGFLAGS)

clean:
	rm -rf daxpy daxpy_c

