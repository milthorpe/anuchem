X10CPP := x10c++
X10FLAG := -O -NO_CHECKS -VERBOSE_CHECKS

ifdef X10RTIMPL # sockets, mpi, pami, etc.
  X10FLAG += -x10rt ${X10RTIMPL}
endif

topology:	Topology.x10
	$(X10CPP) $(X10FLAG) Topology.x10 -o topology

pong:	Pong.x10
	$(X10CPP) $(X10FLAG) Pong.x10 -o pong

clean:
	rm -f topology pong *.cc *.h
