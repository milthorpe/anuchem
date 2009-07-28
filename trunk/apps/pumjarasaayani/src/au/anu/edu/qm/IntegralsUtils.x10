/**
 * IntegralsUtils.x10
 *
 * Some utilities for integral evaluation
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

public class IntegralsUtils { 
    public def this() { 
    } 

    public static def gaussianProductCenter(alpha1:Double, a:Atom, alpha2:Double, b:Atom) : Atom {
        val gamma:Double = alpha1 + alpha2;
        return new Atom(
                         (alpha1 * a.getX() + alpha2 * b.getX()) / gamma,
                         (alpha1 * a.getY() + alpha2 * b.getY()) / gamma,
                         (alpha1 * a.getZ() + alpha2 * b.getZ()) / gamma
                       );
    }

    /**
     * Indexing (i,j,k,l) into long array.
     */
    public static def ijkl2intindex(var i:Int, var j:Int, var k:Int, var l:Int) : Int {
        var temp:Int;

        if (i<j) {
            temp = i;
            i    = j;
            j    = temp;
        } // end if
        if (k<l) {
            temp = k;
            k    = l;
            l    = temp;
        } // end if

        var ij:Int = i * (i+1) / 2+j;
        var kl:Int = k * (k+1) / 2+l;

        if (ij < kl) {
            temp = ij;
            ij   = kl;
            kl   = temp;
        } // end id

        return (ij * (ij+1) / 2+kl);
    }

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

                for(var i:Int=0; i<MAX_ITERATION; i++) {
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

            for(var i:Int=1; i<(MAX_ITERATION+1); i++) {
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
        // TODO: this needs to be changed
        val cof = Array.make[Double]([
                    76.18009172947146,     -86.50532032941677,
                    24.01409824083091,     -1.231739572450155,
                    0.1208650973866179e-2, -0.5395239384953e-5
                  ]);

        var y:Double = x;
        var tmp:Double = x + 5.5;

        tmp -= (x + 0.5) * Math.log(tmp);

        var ser:Double = 1.000000000190015;
        for(var j:Int=0; j<6; j++) {
            y++;
            ser += cof(j) / y;
        } // end for

        return (-tmp + Math.log((2.5066282746310005 * ser) / x));
    }

    // for gammp mathod
    private static SMALL:Double = 0.00000001;
    private static EPS:Double   = 3.0e-7;
    private static FPMIN:Double = 1.0e-30;

    private static MAX_ITERATION:Int = 100;

    /** 
    private static cof:Array[Double] = Array.make[Double]([
        76.18009172947146,     -86.50532032941677,
        24.01409824083091,     -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5
    ]);
    **/
}

