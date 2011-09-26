X10CPP := x10c++
X10CPPFLAGS := -O -x10rt mpi -VERBOSE_CHECKS
X10SOURCEPATH := -sourcepath ../../anu-chem/src
ifdef CPUPROFILE
  X10CPPFLAGS += -post '\# -g -lprofiler'
endif
X10FLAGS := X10_NTHREADS=1

test:	boltzmann
	$(X10FLAGS) ./boltzmann

boltzmann:	Boltzmann.x10
	$(X10CPP) $(X10CPPFLAGS) $(X10SOURCEPATH) Boltzmann.x10 -o boltzmann

clean:
	rm -f boltzmann *.cc *.h
