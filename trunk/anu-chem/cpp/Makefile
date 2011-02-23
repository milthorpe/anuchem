CPP := g++
CFLAGS := -O

run:	directcpp
	./directcpp

directcpp:	ElectrostaticDirect.cc
	$(CPP) $(CFLAGS) ElectrostaticDirect.cc -o directcpp

clean:
	rm -f directcpp
