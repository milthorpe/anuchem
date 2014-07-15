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
 * The range over which to compute the set is distributed across all places.
 * For escaping points outside the set, computes a "fractional iteration count".
 * @see http://linas.org/art-gallery/escape/escape.html
 * @author milthorpe 01/2011
 */
public class MandelbrotDist {
    public static val LIMIT = 2;
    public static val MAX_ITERATIONS = 255;

	public def compute(min : Complex, max : Complex, realPoints : Int) = {
        val gridSpacing = (max.re - min.re) / realPoints;
        val imaginaryPoints = ((max.im - min.im) / (max.re - min.re) * realPoints) as Int;

        val perPlace = realPoints / Place.numPlaces();
        val leftOver = realPoints % Place.numPlaces();
        val firstChunk = leftOver*(perPlace+1);
        val offsets = new Array[Int](Place.numPlaces()+1, (i : Int) => (i < leftOver ? i*(perPlace+1) : firstChunk+(i-leftOver)*perPlace));
        val times = DistArray.make[Double](Dist.makeUnique());
        val result = DistArray.make[Array[Double](2)](Dist.makeUnique());

        for (i in 0..15) {
            Console.OUT.println(offsets);
            finish ateach(place in Dist.makeUnique()) {
                val myResult = new Array[Double](offsets(here.id)..offsets(here.id+1) * 0..(imaginaryPoints-1));
                val start = System.nanoTime();
                for ([gridRe,gridIm] in myResult) {
                    val c = Complex(min.re + gridRe * gridSpacing, min.im + gridIm * gridSpacing);
                    var zn : Complex = c;
                    var i : Int = 0;
                    while (zn.abs() <= LIMIT && i < MAX_ITERATIONS) {
                        zn = zn*zn + c;
                        i++;
                    }
                    if (i == MAX_ITERATIONS) {
                        // series converges at this point, it is part of the Mandelbrot set
                        myResult(gridRe,gridIm) = 0.0;
                    } else {
                        zn = zn*zn + c;
                        zn = zn*zn + c;
                        myResult(gridRe,gridIm) = (i+2) - (Math.log(Math.log(zn.abs())))/ Math.log(2.0);
                    }
                }
                val stop = System.nanoTime();
                times(here.id) = ((stop-start) as Double);
                
                result(here.id) = myResult;
            }

            val increment = Math.max(1.0, realPoints / (Place.numPlaces()*15)) as Int;
            var prevOffset : Int = 0;
            for (place in Place.places()) {
                val thisPlaceTime = at(place) times(here.id);
                Console.OUT.print(thisPlaceTime / 1e6 + " ");
                val tempOffset = offsets(place.id);
                if (place.id > 0) {
                    val share = offsets(place.id) - prevOffset;
                    //Console.OUT.println((place.id-1) + " before, share was " + share);

                    val prevPlaceTime = at(place.prev()) times(here.id);
                    val totalTime = thisPlaceTime + prevPlaceTime;
                    val totalWork = offsets(place.id+1) - prevOffset;
                    val prevShare = ((prevPlaceTime / totalTime) * totalWork as Double) as Int;
                    //Console.OUT.println("previous time " + (prevPlaceTime / totalTime));
                    if ((prevPlaceTime / totalTime) < 0.45) {
                        offsets(place.id) +=  increment;
                    } else if ((prevPlaceTime / totalTime) > 0.55) {
                        offsets(place.id) -=  increment;
                    }
                }
                prevOffset = tempOffset;
            }
            
            Console.OUT.println("");
        }
Console.OUT.println(offsets);
/*
        for (place in Place.places()) at(place) {
            val myResult = result(here.id);
            val pointsHere = myResult.region;
            for ([gridRe] in pointsHere.projection(0)) {
                for ([gridIm] in pointsHere.projection(1)) {
                    Console.OUT.print(myResult(gridRe,gridIm) + " ");
                }
                Console.OUT.println("");
            }
        }
*/
	}

	public static def main(var args: Array[String](1)): void = {
        val min = Complex(-2.0, -1.0);
        val max = Complex(1.0, 1.0);
        var realPoints : Int = 1200;
        if (args.size > 0) {
            realPoints = Int.parseInt(args(0));
        }
		new MandelbrotDist().compute(min, max, realPoints);
	}

}
