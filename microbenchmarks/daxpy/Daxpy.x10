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
public class Daxpy(N:Int) {
    private static ITERS = 1000;

    public def this(N:Int) {
        property(N);
    }

	public def testAll() {
        var start:Long;
        var stop:Long;

        val alpha = 2.5;
        val x = new Rail[Double](N, (i:Long) => i as Double);
        val y = new Rail[Double](N, (i:Long) => i as Double);
        start = System.nanoTime();
        for (iter in 1..ITERS) {
            for (i in 0..(N-1)) {
                x(i) = alpha * x(i) + y(i);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("X10 DAXPY for vectors length %d: %g ms\n", N, ((stop-start) as Double) / 1e6 / ITERS);
	}

	public static def main(args:Rail[String]): void = {
        var size:Int = 100000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Int.parse(args(0));
        }
		new Daxpy(size).testAll();
	}
}
