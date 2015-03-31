/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */

import x10.array.Array_2;
import x10.array.DenseIterationSpace_2;

/**
 * Implements a five-point Laplacian 2D array stencil,
 * similar to the HeatTransfer example from the X10 2.1 tutorial
 * @see http://x10.sourceforge.net/tutorials/x10-2.1/SC_2010/SC10_tut143_X10_Tutorial_final_v3.html
 * @author milthorpe 05/2011
 */
public class FivePointStencil(N:Long) {
    static val EPSILON = 1.0e-5;
    static val TIME_STEPS = 1000;
    val innerRegion:DenseIterationSpace_2;
    var current:Array_2[Double];
    var previous:Array_2[Double];

    public def this(N:Long) {
        property(N);
        innerRegion  = new DenseIterationSpace_2(1, N-2, 1, N-2);
        current = new Array_2[Double](N-1, N-1);
        previous = new Array_2[Double](N-1, N-1);
    }

    /**
     * Performs a five-point stencil on the <code>previous</code> values,
     * storing the result in <code>current</code>,
     * and returns the maximum change in value between previous and current
     * for any point.
     */
    def stencil():Double {
        val temp = previous;
        previous = current;
        current = temp;
        var maxDelta:Double = 0;
        for ([x,y] in innerRegion) {
            current(x,y) = (previous(x, y+1)
                        + previous(x, y-1)
                        + previous(x+1, y)
                        + previous(x-1, y)) / 4.0;
            maxDelta = Math.max(maxDelta, Math.abs(current(x,y) - previous(x,y)));
        }
        return maxDelta;
    }

    /**
     * A single-place parallel version of <code>stencil()</code>.
     */
    def parallelStencil():Double {
        val temp = previous;
        previous = current;
        current = temp;
        var maxDelta:Double = 0;
        finish for (x in 1..(N-2)) async {
            var localMax:Double = 0.0;
            for (y in 1..(N-2)) {
                current(x,y) = (previous(x, y+1)
                            + previous(x, y-1)
                            + previous(x+1, y)
                            + previous(x-1, y)) / 4.0;
                localMax = Math.max(localMax, Math.abs(current(x,y) - previous(x,y)));
            }
            atomic maxDelta = Math.max(localMax, maxDelta);
        }
        return maxDelta;
    }

    /**
     * A single-place parallel version of <code>stencil()</code>,
     * using collecting finish for the max change in temp.
     */
    def parallelCollectingStencil():Double {
        val temp = previous;
        previous = current;
        current = temp;
        val maxDelta = finish (MaxReducer()) {
            for (x in 1..(N-2)) async {
                var localMax:Double = 0.0;
                for (y in 1..(N-2)) {
                    current(x,y) = (previous(x, y+1)
                                + previous(x, y-1)
                                + previous(x+1, y)
                                + previous(x-1, y)) / 4.0;
                    localMax = Math.max(localMax, Math.abs(current(x,y) - previous(x,y)));
                }
                offer localMax;
            }
        };
        return maxDelta;
    }

    static struct MaxReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = Math.max(Math.abs(a),Math.abs(b));
    }

    def initialise() {
        val topEdge = new DenseIterationSpace_2(0, 0, 0, N-1);
        for(p in topEdge) {
            current(p) = 1.0;
            previous(p) = 1.0;
        }
        for (p in innerRegion) {
            current(p) = 0.0;
        }
    }

    def printGrid() {
        Console.OUT.println("grid = ");
        for (i in 0..(N-1)) {
            for (j in 0..(N-1)) {
                Console.OUT.print(current(i,j) + " ");
            }
            Console.OUT.println();
        }
    }

	public def run(): Boolean {
        initialise();
        val start = System.nanoTime();
        for (t in 1..TIME_STEPS) {
            val maxDelta = stencil();
            if (maxDelta < EPSILON) {
                Console.OUT.println("converged after " + t + " time steps");
                break;
            }
        }
        val stop = System.nanoTime();
        //printGrid();
        Console.OUT.printf("sequential five-point stencil avg: %g ms\n", ((stop-start) as Double) / 1e06 / TIME_STEPS);

        initialise();
        val start2 = System.nanoTime();
        for (t in 1..TIME_STEPS) {
            val maxDelta = parallelStencil();
            if (maxDelta < EPSILON) {
                Console.OUT.println("converged after " + t + " time steps");
                break;
            }
        }
        val stop2 = System.nanoTime();
        Console.OUT.printf("parallel five-point stencil avg: %g ms\n", ((stop2-start2) as Double) / 1e06 / TIME_STEPS);

        initialise();
        val start3 = System.nanoTime();
        for (t in 1..TIME_STEPS) {
            val maxDelta = parallelCollectingStencil();
            if (maxDelta < EPSILON) {
                Console.OUT.println("converged after " + t + " time steps");
                break;
            }
        }
        val stop3 = System.nanoTime();
        Console.OUT.printf("parallel five-point stencil with collecting finish avg: %g ms\n", ((stop3-start3) as Double) / 1e06 / TIME_STEPS);

        return true;
	}

	public static def main(args:Rail[String]):void {
        var elementsPerPlace:Long = 512;
        if (args.size > 0) {
            elementsPerPlace = Long.parse(args(0));
        }
		new FivePointStencil(elementsPerPlace).run();
	}
}
