CPP := mpicxx
CFLAGS := -O3 -funroll-loops -g
#LDFLAGS := -lpapi

run:	directcpp
	./directcpp

directcpp:	ElectrostaticDirect.cc
	$(CPP) $(CFLAGS) ElectrostaticDirect.cc $(LDFLAGS) -o directcpp

clean:
	rm -f directcpp
