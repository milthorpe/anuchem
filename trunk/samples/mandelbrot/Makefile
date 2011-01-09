X10CPP := x10c++
XFLAGS := -O -NO_CHECKS -x10rt mpi

plot:	mandelbrot
	./mandelbrot > mandelbrot.dat
	gnuplot < mandelbrot.plot

mandelbrot:	Mandelbrot.x10
	$(X10CPP) $(XFLAGS) Mandelbrot.x10 -o mandelbrot

clean:
	rm -f mandelbrot *.cc *.h mandelbrot.dat mandelbrot.png
