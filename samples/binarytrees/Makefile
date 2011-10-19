X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS
ifdef CPUPROFILE
  X10CPPFLAGS += -post '\# -g -lprofiler'
endif
X10FLAGS := X10_NTHREADS=1

run:	binarytrees
	$(X10FLAGS) ./binarytrees 12

binarytrees:	binarytrees.x10
	$(X10CPP) $(X10CPPFLAGS) binarytrees.x10 -o binarytrees

clean:
	rm -f binarytrees *.cc *.h
