CILKC := cilkc
CILKCFLAGS := -O2 -D_XOPEN_SOURCE=600

run:    parallelLoop
	./parallelLoop

parallelLoop:	parallelLoop.cilk
	 ${CILKC} ${CILKCFLAGS} parallelLoop.cilk -o parallelLoop

clean:
	rm -f parallelLoop
