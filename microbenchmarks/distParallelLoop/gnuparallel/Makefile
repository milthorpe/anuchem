CPP := g++
CPPFLAGS := -O3 -fopenmp -lgomp -std=c++0x

run:    parallelLoop
	./parallelLoop

parallelLoop:	parallelLoop.cc
	${CPP} ${CPPFLAGS} parallelLoop.cc -o parallelLoop

clean:
	rm -f parallelLoop
