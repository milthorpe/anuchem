/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.edu.anu.qm;

/**
 * Some math functions that are not provided by the x10 library
 *
 * @author: V.Ganesh
 */
public final class MathUtil { 
    /**
     * Returns the coefficient of x^j in the expansion of
     * (x+a)^l * (x+b)^m
     */
    public static def binomialPrefactor(j:Int, l:Int, m:Int, a:Double, b:Double): Double {
        var sum:Double = 0.0;
        for(var t:Int = 0; t<(j+1); t++) {
            if(((j-l) <= t) && (t <= m)) {
                sum += binomial(l, j-t) * binomial(m, t)
                      * Math.pow(a, l-j+t) * Math.pow(b, m-t);
            }
        }
        return sum;
    }

    /**
     * Calculates the binomial coefficient
     * (i)      
     * (j) = i! / (i-j)!j!
     * using a simple recurrence relation.  See e.g.
     * http://blog.plover.com/math/choose.html
     */
    public static def binomial(i:Int, j:Int) : Long {
        assert (i > j);
        var n : Long = i;
        var r : Long = 1;
        for (var d:Long=1; d <= j; d++) {
            r *= n--;
            r /= d;
        }
        return r;
    } 

    /**
     * Returns the factorial n!
     */
    public static def factorial(var n:Int) : Long {
        var value:Long = 1;
        
        while(n > 1) {
            value = value * n;
            n--;
        } // end while
        
        return value;
    }

    /**
     * Returns the double factorial n!! =
     * 1,          if n==0 || n==1
     * n((n-2)!),  if n > 1
     */
    public static def factorial2(var n:Int) : Long {
        var value:Long = 1;
        
        while(n > 0) {
            value = value * n;
            n-=2;
        } // end while
        
        return value;
    }

    /**
     * Calculates the ratio
     * i! / (i-2*j)!j!
     * using a simple recurrence relation.  See e.g.
     * http://blog.plover.com/math/choose.html
     */
    public static def factorialRatioSquared(a:Int, b:Int) : Double {
        assert (a > b);
        var n : Long = a;
        var r : Long = 1;
        for (var d:Long=1; d <= b; d++) {
            r *= n--;
            r /= d;
        }
        for (var i:Int=(a-b); i>(a-2*b); i--) {
            r *= i;
        }
        return r;
    }

    /**
     * @return minus one to the power (a)
     */
    public static final def minusOnePow(a: Int) {
        return (a % 2 == 0) ? 1 : -1;
    }
}

