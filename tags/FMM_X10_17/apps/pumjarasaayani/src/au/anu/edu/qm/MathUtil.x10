/**
 * MathUtil.x10
 *
 * Some math functions that are not provided by the x10 library
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

public final class MathUtil { 
    public def this() { 
    } 

    public static def binomialPrefactor(s:Int, ia:Int, ib:Int, xpa:Double, xpb:Double): Double {
        var sum:Double = 0.0;

        for(var t:Int = 0; t<(s+1); t++) {
            if(((s-ia) <= t) && (t <= ib)) {
                sum += binomial(ia, s-t) * binomial(ib, t)
                      * Math.pow(xpa, ia-s+t) * Math.pow(xpb, ib-t);
            } // end if
        } // end for

        return sum;
    }
   
    public static def binomial(i:Int, j:Int) : Double {
        return (factorial(i) / factorial(j) / factorial(i - j));
    } 

    public static def factorial(var n:Int) : Long {
        var value:Long = 1;
        
        while(n > 1) {
            value = value * n;
            n--;
        } // end while
        
        return value;
    }

    public static def factorial2(var n:Int) : Long {
        var value:Long = 1;
        
        while(n > 0) {
            value = value * n;
            n-=2;
        } // end while
        
        return value;
    }

    public static def factorialRatioSquared(a:Int, b:Int) : Double {
        return factorial(a) / factorial(b) / factorial(a-2*b);
    }
}

