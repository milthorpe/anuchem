CPP := mpicxx
CFLAGS := -O2 -g
LDFLAGS := -lpapi

run:	directcpp
	./directcpp

directcpp:	ElectrostaticDirect.cc
	$(CPP) $(CFLAGS) ElectrostaticDirect.cc $(LDFLAGS) -o directcpp

clean:
	rm -f directcpp
