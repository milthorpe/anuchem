X10CPP := x10c++
X10RT := sockets
X10CPPFLAGS := -O -NO_CHECKS -x10rt $(X10RT) -VERBOSE_CHECKS -x10lib ../../xla/xla_$(X10RT).properties -x10lib ../../anu-chem/anuchem_$(X10RT).properties
ifdef CPUPROFILE
  X10CPPFLAGS += -gpt -DEBUG -post '\# -L/scratch/perftools/lib'
endif

test:	cutoff
	./cutoff 100000

cutoff:	ElectrostaticCutoff.x10
	$(X10CPP) $(X10CPPFLAGS) ElectrostaticCutoff.x10 -o cutoff

clean:
	rm -f cutoff *.cc *.h
