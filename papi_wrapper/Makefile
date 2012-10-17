CCFLAGS := -O3 -g
LDFLAGS := -lpapi

measure_sqrt:   measure_sqrt.cc
	g++ $(CCFLAGS) measure_sqrt.cc $(LDFLAGS) -o measure_sqrt

measure_div:   measure_div.cc
	g++ $(CCFLAGS) measure_div.cc $(LDFLAGS) -o measure_div

clean:
	rm -f measure_sqrt measure_div
