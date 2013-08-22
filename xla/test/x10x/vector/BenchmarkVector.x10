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

package x10x.vector;

import x10.util.Random;
import harness.x10Test;

/**
 * Benchmarks XLA Vector methods
 * @author milthorpe
 */
public class BenchmarkVector extends x10Test {
    val N:Long;
    val rand = new Random(47);

    public def this(N:Long) {
        super();
        this.N = N;
    }

    public def run(): boolean {
        Console.OUT.println("vector size " + N);
        val a = randomVector();

        var start : Long = System.nanoTime();
        val x = a.magnitude();
        var stop : Long = System.nanoTime();
        Console.OUT.printf("magnitude: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val y = a.maxNorm();
        stop = System.nanoTime();
        Console.OUT.printf("maxNorm: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val a2 = a.normalize();
        stop = System.nanoTime();
        Console.OUT.printf("normalize: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val a3 = a.negate();
        stop = System.nanoTime();
        Console.OUT.printf("negate: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val a4 = a.mul(2.0);
        stop = System.nanoTime();
        Console.OUT.printf("mul: %g ms\n", ((stop-start) as Double) / 1e6);

        val b = randomVector();
        start = System.nanoTime();
        val a5 = a.add(b);
        stop = System.nanoTime();
        Console.OUT.printf("add: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val a6 = a.sub(b);
        stop = System.nanoTime();
        Console.OUT.printf("sub: %g ms\n", ((stop-start) as Double) / 1e6);

        start = System.nanoTime();
        val z = a.dot(b);
        stop = System.nanoTime();
        Console.OUT.printf("dot: %g ms\n", ((stop-start) as Double) / 1e6);

        return true;
    }

    private def randomVector() {
        val a = new Vector(N);
        for (i in 0..(a.vec.size-1)) {
            a(i) = rand.nextDouble();
        }
        return a;
    }

    public static def main(args:Rail[String]) {
        var N:Long = 10;
        if (args.size > 0) {
            N = Long.parseLong(args(0));
        }
        new BenchmarkVector(N).execute();
    }
}
