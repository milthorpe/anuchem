/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */
package x10x.matrix;

import x10.util.Random;
import harness.x10Test;

/**
 * Benchmarks XLA Matrix methods
 * @author milthorpe
 */
public class BenchmarkMatrix extends x10Test {
    val N:Long;
    val rand = new Random(47);

    public def this(N:Long) {
        super();
        this.N = N;
    }

    public def run(): boolean {
        Console.OUT.println("square matrix size " + N);
        val a = randomMatrix();
        val start = System.nanoTime();
        val x = a.sumOffDiagonal();
        val stop = System.nanoTime();
        Console.OUT.printf("sumOffDiagonal: %g ms\n", ((stop-start) as Double) / 1e6);

        val start2 = System.nanoTime();
        val a2 = a.mul(2.0);
        val stop2 = System.nanoTime();
        Console.OUT.printf("scale: %g ms\n", ((stop2-start2) as Double) / 1e6);

        val start3 = System.nanoTime();
        val c = a.transpose();
        val stop3 = System.nanoTime();
        Console.OUT.printf("transpose: %g ms\n", ((stop3-start3) as Double) / 1e6);

        val b = randomMatrix();
        val start4 = System.nanoTime();
        val d = a.mul(b);
        val stop4 = System.nanoTime();
        Console.OUT.printf("gemm: %g ms\n", ((stop4-start4) as Double) / 1e6);

        val start5 = System.nanoTime();
        val e = a.symmetricOrthogonalization();
        val stop5 = System.nanoTime();
        Console.OUT.printf("orthog: %g ms\n", ((stop5-start5) as Double) / 1e6);

        return true;
    }

    private def randomMatrix() {
        val a = new Matrix(N);
        for ([i,j] in a.mat.indices()) {
            a(i,j) = rand.nextDouble();
        }
        return a;
    }

    public static def main(args:Rail[String]) {
        var N:Long = 70;
        if (args.size > 0) {
            N = Long.parseLong(args(0));
        }
        new BenchmarkMatrix(N).execute();
    }
}
