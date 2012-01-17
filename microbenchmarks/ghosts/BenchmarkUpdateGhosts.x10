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
        val r = 0..(arrayDim*facI-1) * 0..(arrayDim*facJ-1) * 0..(arrayDim-1);
        Console.OUT.println("r =  " + r);
        val d = Dist.makeBlockBlock(r, 0, 1);
        //Console.OUT.println("d = " + d);

        val a = DistArray.make[Double](d, 2, false);

        val shiftTime = finish(MaxReducer()) {
            ateach(p in Dist.makeUnique(d.places())) {
                val startHere = System.nanoTime();
                for ([t] in 1..ITERS) {
                    a.sendGhosts();
                    a.waitOnGhosts();
                }
                val stopHere = System.nanoTime();
                offer (stopHere-startHere);
            }
        };
        Console.OUT.printf("updateGhosts shift avg: %g ms\n", ((shiftTime) as Double) / (1e06 * ITERS));

        val b = DistArray.make[Double](d, 2, false, false);
        val putTime = finish(MaxReducer()) {
            ateach(p in Dist.makeUnique(d.places())) {
                val startHere = System.nanoTime();
                for ([t] in 1..ITERS) {
                    b.sendGhosts();
                    b.waitOnGhosts();
                }
                val stopHere = System.nanoTime();
                offer (stopHere-startHere);
            }
        };
        Console.OUT.printf("updateGhosts put avg: %g ms\n", ((putTime) as Double) / (1e06 * ITERS));

        return true;
    }

    static struct MaxReducer implements Reducible[Long] {
        public def zero() = -1L;
        public operator this(a:Long, b:Long) = Math.max(a,b);
    }

    public static def main(var args: Rail[String]): void = {
        var arrayDim : Int = 100;
        if (args.size > 0) {
            arrayDim = Int.parse(args(0));
        }
        new BenchmarkUpdateGhosts(arrayDim).run();
    }

}
