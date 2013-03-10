/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2013.
 */

package org.netlib.fdlibm;

/** 
 * This class prints values of erf(x) and erfc(x) in the range [0.0,5.0] using
 * the X10 Math functions, and compares with results from the Netlib FDLIBM
 * implementation.  The error (difference between Netlib Erf and X10/C++ Math)
 * is printed for each value.
 */
public class TestErf {
    static POINTS = 100;
    public static def main(args:Rail[String]) {
        Console.OUT.println("          Math.erf       err(Erf.erf)          Math.erfc      err(Erf.erfc)");
        for(i in 0..POINTS) {
            val x = i * 5.0 / POINTS;
            Console.OUT.printf("%18e %18e %18e %18e\n", Math.erf(x), Math.erf(x)-Erf.erf(x), Math.erfc(x), Math.erfc(x)-Erf.erfc(x));
        }
        at (here.next()) { }
/*
        var start:Long = System.nanoTime();
        for (iter in 1..100) {
            for (i in 0..POINTS) {
                val x = i * 5.0 / POINTS;
                val e = Math.erfc(x);
            }
        }
        var stop:Long = System.nanoTime();
        Console.OUT.println(stop-start);

        start = System.nanoTime();
        for (iter in 1..100) {
            for (i in 0..POINTS) {
                val x = i * 5.0 / POINTS;
                val e = Erf.erfc(x);
            }
        }
        stop = System.nanoTime();
        Console.OUT.println(stop-start);
*/
    }
}
