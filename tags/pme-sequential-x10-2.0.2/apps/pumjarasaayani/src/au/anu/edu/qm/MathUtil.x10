/**
 * MathUtil.x10
 *
 * Some math functions that are not provided by the x10 library
 *
 * @author: V.Ganesh
 */
package au.anu.edu.qm;

public final class MathUtil { 
    /**
     * Returns the coefficient of x^j in the expansion of
     * (x+a)^l * (x+b)^m
     */
    public static safe def binomialPrefactor(j:Int, l:Int, m:Int, a:Double, b:Double): Double {
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
     * Returns the binomial coefficient
     * (i)      
     * (j) = i! / (i-j)!j!
     */
    public static safe def binomial(i:Int, j:Int) : Double {
        return (factorial(i) / factorial(j) / factorial(i - j));
    } 

    /**
     * Returns the factorial n!
     */
    public static safe def factorial(var n:Int) : Long {
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
    public static safe def factorial2(var n:Int) : Long {
        var value:Long = 1;
        
        while(n > 0) {
            value = value * n;
            n-=2;
        } // end while
        
        return value;
    }

    public static safe def factorialRatioSquared(a:Int, b:Int) : Double {
        return factorial(a) / factorial(b) / factorial(a-2*b);
    }
}

