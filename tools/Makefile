X10CPP := x10c++
XFLAGS := -O -NO_CHECKS

topology:	Topology.x10
	$(X10CPP) $(XFLAGS) Topology.x10 -o topology

pong:	Pong.x10
	$(X10CPP) $(XFLAGS) Pong.x10 -o pong

clean:
	rm -f topology pong *.cc *.h
