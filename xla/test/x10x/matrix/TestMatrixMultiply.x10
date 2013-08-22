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
 * Tests XLA Matrix class multiplication method
 * @author milthorpe
 */
public class TestMatrixMultiply extends x10Test {
    val N:Long;
    val rand = new Random(47);

    public def this(N:Long) {
        super();
        this.N = N;
    }

    public def run(): boolean {
        val a = randomMatrix();
        Console.OUT.println(a);
        val b = randomMatrix();
        Console.OUT.println(b);

        val c = a.mul(b);

        Console.OUT.println(c);

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
        var N:Long = 25;
        if (args.size > 0) {
            N = Long.parseLong(args(0));
        }
        new TestMatrixMultiply(N).execute();
    }
}
