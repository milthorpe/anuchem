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
 * Tests performance of DistArray.updateGhosts
 * N.B. ghost updates not yet integrated to X10 SVN
 * @author milthorpe 09/2011
 */
public class BenchmarkUpdateGhosts(arrayDim : Int) {
    public static ITERS = 10000;

    public def this(elementsPerPlace : Int) {
        property(elementsPerPlace);
    }

    public def run(): Boolean = {
        val facI = Math.sqrt(Place.MAX_PLACES) as Int;
        val facJ = Place.MAX_PLACES / facI;
        val r = 0..513 * 0..513 * 0..9;
        Console.OUT.println("r =  " + r);
        val d = Dist.makeBlockBlock(r, 0, 1);
        //Console.OUT.println("d = " + d);

        val a = DistArray.make[Double](d, 2);
/*
        var start:Long = System.nanoTime();
        for ([t] in 1..ITERS) {
            a.updateGhosts();
        }
        var stop:Long = System.nanoTime();

        Console.OUT.printf("updateGhosts put avg: %g ms\n", ((stop-start) as Double) / (1e06 * ITERS));
/*
        start = System.nanoTime();
        for ([t] in 1..ITERS) {
            a.updateGhostsShift();
        }
        stop = System.nanoTime();

        Console.OUT.printf("updateGhosts shift avg: %g ms\n", ((stop-start) as Double) / (1e06 * ITERS));
*/
		val maxTime = new Accumulator[Long](Reducible.MaxReducer[Long](0L));
        finish ateach(p in Dist.makeUnique(d.places())) {
            val startHere = System.nanoTime();
            for ([t] in 1..ITERS) {
                a.sendGhosts();
                a.waitOnGhosts();
            }
            val stopHere = System.nanoTime();
            maxTime <- (stopHere-startHere);
        }
        Console.OUT.printf("updateGhosts no sync avg: %g ms\n", ((maxTime()) as Double) / (1e06 * ITERS));

        return true;
    }

    public static def main(var args: Rail[String]): void = {
        var arrayDim : Int = 100;
        if (args.size > 0) {
            arrayDim = Int.parse(args(0));
        }
        new BenchmarkUpdateGhosts(arrayDim).run();
    }

}
