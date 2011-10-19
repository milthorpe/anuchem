X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS
ifdef CPUPROFILE
  X10CPPFLAGS += -post '\# -g -lprofiler'
endif
X10FLAGS := X10_NTHREADS=1

plot:	heat
	$(X10FLAGS) ./heat 512

heat:	FivePointStencil.x10
	$(X10CPP) $(X10CPPFLAGS) FivePointStencil.x10 -o heat

clean:
	rm -f heat *.cc *.h
