X10CPP := x10c++
X10RT := sockets
X10CPPFLAGS := -O -DEBUG -NO_CHECKS -x10rt $(X10RT) -VERBOSE_CHECKS -x10lib ../../xla/xla.properties -x10lib ../../anu-chem/anuchem.properties
ifdef CPUPROFILE
  X10CPPFLAGS += -gpt -DEBUG -post '\# -L/scratch/perftools/lib'
endif

test:	cutoff
	./cutoff 100000 10.0

cutoff:	ElectrostaticCutoff.x10
	$(X10CPP) $(X10CPPFLAGS) ElectrostaticCutoff.x10 -o cutoff

cutoffDist:	ElectrostaticCutoffDist.x10
	$(X10CPP) $(X10CPPFLAGS) ElectrostaticCutoffDist.x10 -o cutoffDist

clean:
	rm -f cutoff cutoffDist *.cc *.h
