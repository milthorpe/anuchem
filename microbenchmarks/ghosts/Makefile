X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -x10rt mpi -VERBOSE_CHECKS
ifdef CPUPROFILE
  X10CPPFLAGS += -gpt -DEBUG -post '\# -L/scratch/perftools/lib'
endif
X10FLAGS := X10_NTHREADS=1 X10_STATIC_THREADS=true

run:    updateGhosts
	${X10FLAGS} ./updateGhosts 1025

updateGhosts:	BenchmarkUpdateGhosts.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkUpdateGhosts.x10 -o updateGhosts

clean:
	rm -rf updateGhosts *.cc *.h *.dat *.ps

