CPP := mpicxx
CFLAGS := -O2 -funroll-loops

run:	directcpp
	./directcpp

directcpp:	ElectrostaticDirect.cc
	$(CPP) $(CFLAGS) ElectrostaticDirect.cc -o directcpp

clean:
	rm -f directcpp
