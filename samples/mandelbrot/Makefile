X10CPP := x10c++
XFLAGS := -O -NO_CHECKS

plot:	mandelbrot
	./mandelbrot > mandelbrot.dat
	gnuplot < mandelbrot.plot

mandelbrot:	Mandelbrot.x10
	$(X10CPP) $(XFLAGS) Mandelbrot.x10 -o mandelbrot

mandelbrot_dist:	MandelbrotDist.x10
	$(X10CPP) $(XFLAGS) -x10rt mpi Mandelbrot.x10 -o mandelbrot

clean:
	rm -f mandelbrot *.cc *.h mandelbrot.dat mandelbrot.png
