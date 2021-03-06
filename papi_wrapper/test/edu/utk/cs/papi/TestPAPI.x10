/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package edu.utk.cs.papi;

/**
 * Tests wrapper for PAPI toolkit
 * Counts FLOPs and memory ops for a dot product
 * @see http://icl.cs.utk.edu/papi
 * @author milthorpe
 */
public class TestPAPI {
    public static def main(args:Rail[String]) {
        val N = 1000000;
        val a = new Array[Double](N, 1.00001);

        val papi = new PAPI();
        papi.countFlops();

        papi.start();
        var x:Double=1.0;
        for (i in 0..(N-1)) {
            x *= a(i);
        }
        papi.stop();
        papi.start();
        for (i in 0..(N-1)) {
            x *= a(i);
        }
        papi.stop();
        papi.printFlops();

        Console.OUT.println("dot product = " + x);

        papi.countMemoryOps();

        papi.start();
        x=1.0;
        for (i in 0..(N-1)) {
            x *= a(i);
        }
        papi.stop();
        papi.printMemoryOps();

        Console.OUT.println("dot product = " + x);
    }
}

