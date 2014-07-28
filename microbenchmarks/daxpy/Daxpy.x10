/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2013.
 */
import x10.compiler.Inline;

/**
 * Benchmarks simple DAXPY operation
 * @author milthorpe 07/2013
 */
public class Daxpy(N:Long) {
    private static ITERS = 1000;

    public def this(N:Long) {
        property(N);
    }

    private def daxpy(alpha:Double, x:Rail[Double], y:Rail[Double]) {
        val start = 0;
        val end = x.size-1;
/*
        for (i in start..end) {
            x(i) = alpha * x(i) + y(i);
        }
*/
        val body = (i:Long) => {
            x(i) = alpha * x(i) + y(i);
        };
/*
        finish for (i in start..end) async {
            body(i);
        }
*/

        val numElem = end - start + 1;
        val blockSize = numElem / Runtime.NTHREADS;
        val leftOver = numElem % Runtime.NTHREADS;
        finish for (var t:Long=Runtime.NTHREADS-1; t>=0; t--) {
            val tStart = start + t <= leftOver ? t*(blockSize+1) : t*blockSize + leftOver;
            val tEnd = tStart + ((t < leftOver) ? (blockSize+1) : blockSize);
            async {
                for (i in tStart..tEnd) {
                    body(i);
                }
            }
        }

/*
        finish RecursiveBisection1D(start, end+1, 16).execute(body);
*/
    }

	public def testAll() {
        var start:Long;
        var stop:Long;

        val alpha = 2.5;
        val x = new Rail[Double](N, (i:Long) => i as Double);
        val y = new Rail[Double](N, (i:Long) => i as Double);
        start = System.nanoTime();
        for (iter in 1..ITERS) {
            daxpy(alpha, x, y);
        }
        stop = System.nanoTime();
        Console.OUT.printf("X10 DAXPY for vectors length %d: %g ms\n", N, ((stop-start) as Double) / 1e6 / ITERS);
	}

	public static def main(args:Rail[String]): void = {
        var size:Long = 100000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Long.parse(args(0));
        }
		new Daxpy(size).testAll();
	}

    private static struct RecursiveBisection1D(start:Long, end:Long, grainSize:Long) {
        public def this(start:Long, end:Long) {
            val grainSize = (end-start) / (Runtime.NTHREADS*8);
            property(start, end, grainSize);
        }

        public def this(start:Long, end:Long, grainSize:Long) {
            property(start, end, grainSize);
        }

        public def execute(body:(idx:Long)=> void) {
            if ((end-start) > grainSize) {
                val secondHalf=RecursiveBisection1D((start+end)/2L, end, grainSize);
                async secondHalf.execute(body);
                val firstHalf=RecursiveBisection1D(start, (start+end)/2L, grainSize);
                firstHalf.execute(body);
            } else {
                for (i in start..(end-1)) {
                    body(i);
                }
            }
        }
    }
}
