/*
 *  This file is part of the X10 project (http://x10-lang.org).
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2011.
 */

/**
 * Computes the Mandelbrot set for the given range of complex numbers.
 * For escaping points outside the set, computes a "fractional iteration count".
 * @see http://linas.org/art-gallery/escape/escape.html
 * @author milthorpe 01/2011
 */
public class Mandelbrot {
    public static val LIMIT = 2;
    public static val MAX_ITERATIONS = 255;

	public def compute(min : Complex, max : Complex, realPoints : Int) = {
        val gridSpacing = (max.re - min.re) / realPoints;
        val imaginaryPoints = ((max.im - min.im) / (max.re - min.re) * realPoints) as Int;
        val result = DistArray.make[Double](Dist.makeBlock(0..realPoints * 0..imaginaryPoints));
        finish ateach (place in Dist.makeUnique()) {
            val start = System.nanoTime();
            for ([gridRe,gridIm] in result.dist(here)) {
                val c = Complex(min.re + gridRe * gridSpacing, min.im + gridIm * gridSpacing);
                var zn : Complex = c;
                var i : Int = 0;
                while (zn.abs() <= LIMIT && i < MAX_ITERATIONS) {
                    zn = zn*zn + c;
                    i++;
                }
                if (i == MAX_ITERATIONS) {
                    // series converges at this point, it is part of the Mandelbrot set
                    result(gridRe,gridIm) = 0.0;
                } else {
                    zn = zn*zn + c;
                    zn = zn*zn + c;
                    result(gridRe,gridIm) = (i+2) - (Math.log(Math.log(zn.abs())))/ Math.log(2.0);
                }
            }
            val stop = System.nanoTime();
            Console.OUT.printf("# time at " + here + " : %g ms\n", ((stop-start) as Double) / 1e6);
        }
        for ([gridRe] in result.region.projection(0)) {
            for ([gridIm] in result.region.projection(1)) {
                at (result.dist(gridRe,gridIm)) {
                    Console.OUT.print(result(gridRe,gridIm) + " ");
                }
            }
            Console.OUT.println("");
        }
	}

	public static def main(var args: Array[String](1)): void = {
        val min = Complex(-2.0, -1.0);
        val max = Complex(1.0, 1.0);
        var realPoints : Int = 640;
        if (args.size > 0) {
            realPoints = Int.parseInt(args(0));
        }
		new Mandelbrot().compute(min, max, realPoints);
	}

}
