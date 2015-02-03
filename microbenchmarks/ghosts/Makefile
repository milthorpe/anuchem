X10CPP := x10c++
X10CPPFLAGS := -O -NO_CHECKS -VERBOSE_CHECKS
ifdef PROFILE
  X10CPPFLAGS += -gpt
endif
X10FLAGS := X10_NPLACES=4 X10_NTHREADS=1

run:    updateGhosts
	${X10FLAGS} ./updateGhosts 128

updateGhosts:	BenchmarkUpdateGhosts.x10
	$(X10CPP) $(X10CPPFLAGS) BenchmarkUpdateGhosts.x10 -o updateGhosts

clean:
	rm -rf updateGhosts *.cc *.h *.dat *.ps

