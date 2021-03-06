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
 * Some utilities for integral evaluation
 *
 * @author: V.Ganesh
 */
public class IntegralsUtils { 
    /**
     * Incomplete gamma function
     */
    public static def computeFGamma(m:Int, var x:Double) : Double {
        x = Math.max(Math.abs(x), SMALL);

        return (0.5 * Math.pow(x, -m - 0.5)
                     * gammaIncomplete(m + 0.5, x));
    }

    /**
     * Incomple gamma function gamma() computed from Numerical Recipes
     * routine gammp.
     */
    public static def gammaIncomplete(a:Double, x:Double) : Double {
        var gammap:Double = 0.0;
        var gln:Double    = 0.0;

        gln = gammln(a);

        if(x < (a + 1.0)) {
            // Series representation of Gamma. NumRec sect 6.1.
            if(x != 0.0) {
                var ap:Double = a;
                var sum:Double;
                var delta:Double = sum = 1.0 / a;

                for(var i:Int=0n; i<MAX_ITERATION; i++) {
                    ap++;
                    delta *= x / ap;
                    sum += delta;
                    if(Math.abs(delta) < Math.abs(sum) * EPS) break;
                } // end for

                gammap = sum * Math.exp(-x + a*Math.log(x) - gln);
            } else {
                gammap = 0.0;
            } // end if
        } else {
            // Continued fraction representation of Gamma. NumRec sect 6.1
            var b:Double = (x + 1.0) - a;
            var c:Double = 1.0 / FPMIN;
            var d:Double = 1.0 / b;
            var h:Double = d;
            var an:Double, delta:Double;

            for(var i:Int=1n; i<(MAX_ITERATION+1n); i++) {
                an = -i * (i-a);
                b += 2.0;
                d = an * d + b;

                if(Math.abs(d) < FPMIN) d = FPMIN;

                c = b + an / c;

                if(Math.abs(c) < FPMIN) c = FPMIN;

                d = 1.0 / d;
                delta = d * c;
                h *= delta;

                if(Math.abs(delta - 1.0) < EPS) break;
            } // end for

            gammap = 1.0 - (Math.exp(-x + a*Math.log(x) - gln) * h);
        } // end if

        return (Math.exp(gln) * gammap);
    }

    /**
     * Numerical recipes, section 6.1
     */
    public static def gammln(x:Double) : Double {
        var y:Double = x;
        var tmp:Double = x + 5.5;

        tmp -= (x + 0.5) * Math.log(tmp);

        var ser:Double = 1.000000000190015;

            y++;
            ser += 76.18009172947146 / y;
            y++;
            ser += -86.50532032941677 / y;
            y++;
            ser += 24.01409824083091 / y;
            y++;
            ser += -1.231739572450155 / y;
            y++;
            ser += 0.1208650973866179e-2 / y;
            y++;
            ser += -0.5395239384953e-5 / y;

        return (-tmp + Math.log((2.5066282746310005 * ser) / x));
    }

    // for gammp mathod
    private static SMALL:Double = 0.00000001;
    private static EPS:Double   = 3.0e-7;
    private static FPMIN:Double = 1.0e-30;

    private static MAX_ITERATION:Int = 100n;
}

