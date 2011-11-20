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
package au.edu.anu.util;

/** A heavyweight timer for use in gathering performance statistics. */
final public class StatisticalTimer {
    public val current:Rail[Long];
    public val count:Rail[Long];
    public val total:Rail[Long];
    public val min:Rail[Long];
    public val max:Rail[Long];
    public val sumOfSquares:Rail[Double];

    public def this(n:Int) {
        current = new Rail[Long](n);
        count = new Rail[Long](n);
        total = new Rail[Long](n);
        min = new Rail[Long](n, Long.MAX_VALUE);
        max = new Rail[Long](n, Long.MIN_VALUE);
        sumOfSquares = new Rail[Double](n);
    }

    public def start(id:Int) { current(id) = -System.nanoTime(); }
    public def clear(id:Int) { 
        total(id) = 0; count(id) = 0;
        min(id) = Long.MAX_VALUE; max(id) = Long.MIN_VALUE;
        sumOfSquares(id) = 0.0;
    }

    public def stop(id:Int) {
        current(id) += System.nanoTime();
        val elapsed = current(id);
        count(id)++;
        total(id) += elapsed;
        min(id) = Math.min(min(id), elapsed);
        max(id) = Math.max(max(id), elapsed);
        val e = elapsed as Double;
        sumOfSquares(id) += e*e;
    }

    public def last(id:Int) = current(id);
    public def mean(id:Int) = total(id) / count(id);
    public def stdDev(id:Int):Double {
        val s0 = count(id) as Double;
        val s1 = total(id) as Double;
        val s2 = sumOfSquares(id);
        return Math.sqrt(s0*s2 - s1*s1) / s0;
    }
}

