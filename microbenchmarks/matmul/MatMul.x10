/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */

//import x10.array.DenseIterationSpace_2;

/**
 * Benchmarks naive matrix multiplication
 * @author milthorpe 07/2013
 */
public class MatMul(N:Long) {
    private static ITERS = 100;

    public def this(N:Long) {
        property(N);
    }

    private def matmul(a:Rail[Double], b:Rail[Double], c:Rail[Double]) {
/* original
        val indices = new DenseIterationSpace_2(0, 0, (N-1), (N-1));
        for ([i,j] in indices) {
            var x:Double = 0.0;
            for (k in 0..(N-1)) {
                x += a(i+k*N) * b(k+j*N);
            }
            c(i+j*N) = x;
        }
*/
        val body = (min_i:Long, max_i:Long, min_j:Long, max_j:Long) => {
            for (j in min_j..max_j) {
                for (i in min_i..max_i) {
                    var x:Double = 0.0;
                    for (k in 0..(N-1)) {
                        x += a(i+k*N) * b(k+j*N);
                    }
                    c(i+j*N) = x;
                }
            }
        };
/*
        body(0, N-1, 0, N-1);
*/
/*
        finish for (i in 0..(N-1)) {
            for (j in 0..(N-1)) {
                async body(i, i, j, j);
            }
        }
*/
/*
        val numElem1 = N;
        val blockSize1 = numElem1 / Runtime.NTHREADS;
        val leftOver1 = numElem1 % Runtime.NTHREADS;
        val numElem2 = N;
        val blockSize2 = numElem2 / Runtime.NTHREADS;
        val leftOver2 = numElem2 % Runtime.NTHREADS;
        finish {
            for (t1 in 0..(Runtime.NTHREADS-1)) {
                for (t2 in 0..(Runtime.NTHREADS-1))  async {
                    val tMin_i1 = 0 + t1 <= leftOver1 ? t1*(blockSize1+1) : t1*blockSize1 + leftOver1;
                    val tMax_i1 = tMin_i1 + ((t1 < leftOver1) ? (blockSize1+1) : blockSize1) - 1;
                    val tMin_i2 = 0 + t2 <= leftOver2 ? t2*(blockSize2+1) : t2*blockSize2 + leftOver2;
                    val tMax_i2 = tMin_i2 + ((t2 < leftOver2) ? (blockSize2+1) : blockSize2) - 1;
                    body(tMin_i1, tMax_i1, tMin_i2, tMax_i2);
                }
            }
        }
*/

        finish RecursiveBisection2D(0, N, 0, N).execute(body);

    }

	public def testAll() {
        var start:Long;
        var stop:Long;

        val a = new Rail[Double](N*N, (i:Long) => i as Double);
        val b = new Rail[Double](N*N, (i:Long) => i as Double);
        val c = new Rail[Double](N*N, (i:Long) => i as Double);
        start = System.nanoTime();
        for (iter in 1..ITERS) {
            matmul(a, b, c);
        }
        stop = System.nanoTime();
        Console.OUT.printf("X10 matrix multiplication size %d: %g ms\n", N, ((stop-start) as Double) / 1e6 / ITERS);
	}

	public static def main(args:Rail[String]): void = {
        var size:Long = 500;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Long.parse(args(0));
        }
		new MatMul(size).testAll();
	}

    private static struct RecursiveBisection2D(s1:Long, e1:Long, s2:Long, e2:Long, g1:Long, g2:Long) {
        public def this(s1:Long, e1:Long, s2:Long, e2:Long) {
            val g1 = (e1-s1) / Runtime.NTHREADS;
            val g2 = (e2-s2) / Runtime.NTHREADS;
            property(s1, e1, s2, e2, g1, g2);
        }

        public def this(s1:Long, e1:Long, s2:Long, e2:Long, g1:Long, g2:Long) {
            property(s1, e1, s2, e2, g1, g2);
        }

        public def execute(body:(min_i1:Long, max_i1:Long, min_i2:Long, max_i2:Long)=> void) {
            if ((e1-s1) > g1 && ((e1-s1) >= (e2-s2) || (e2-s2) <= g2)) {
                val secondHalf=RecursiveBisection2D((s1+e1)/2L, e1, s2, e2, g1, g2);
                async secondHalf.execute(body);
                val firstHalf=RecursiveBisection2D(s1, (s1+e1)/2L, s2, e2, g1, g2);
                firstHalf.execute(body);
            } else if ((e2-s2) > g2) {
                val secondHalf=RecursiveBisection2D(s1, e1, (s2+e2)/2L, e2, g1, g2);
                async secondHalf.execute(body);
                val firstHalf=RecursiveBisection2D(s1, e1, s2, (s2+e2)/2L, g1, g2);
                firstHalf.execute(body);
            } else {
                body(s1, e1-1, s2, e2-1);
            }
        }
    }
}
