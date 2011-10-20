/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 * (C) Copyright Josh Milthorpe 2011.
 */
package au.edu.anu.qm;

import x10x.vector.Point3d;
import x10.compiler.Inline;

/**
 * Represents a primitive gaussian function
 *
 * @author: V.Ganesh
 */
public final class PrimitiveGaussian {
    public val origin:Point3d;
    public val power:Power;
    public val exponent:Double;
    public var coefficient:Double;
    public val normalization:Double;

    public def this(origin:Point3d, power:Power, exponent:Double, coefficient:Double, normalize:Boolean) { 
        this.origin = origin;
        this.power = power;
        this.exponent = exponent;
        this.coefficient = coefficient;
        if (normalize) {
            val l = power.l; 
            val m = power.m;
            val n = power.n;       

            this.normalization = Math.sqrt(Math.pow(2, 2 * (l + m + n) + 1.5) *
                                       Math.pow(exponent, l + m + n + 1.5) /
                                       MathUtil.factorial2(2 * l - 1) / 
                                       MathUtil.factorial2(2 * m - 1) /
                                       MathUtil.factorial2(2 * n - 1) /
                                       PI_RAISE_TO_1DOT5);
        } else {
            this.normalization = 0.0; // TODO - should be 1.0?
        }
    } 

    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    public def setCoefficient(coef:Double) { coefficient = coef; }

    public def overlap(pg:PrimitiveGaussian) : Double {
        return normalization * pg.normalization * ovrlp(pg);
    }

    private def ovrlp(pg:PrimitiveGaussian) : Double {
        val radiusABSquared = origin.distanceSquared(pg.origin);
        val prod = mul(pg);

        val wx = overlap1D(power.l, pg.power.l,
                            prod.origin.i - origin.i,
                            prod.origin.i - pg.origin.i, prod.exponent);
        val wy = overlap1D(power.m, pg.power.m,
                            prod.origin.j - origin.j,
                            prod.origin.j - pg.origin.j, prod.exponent);
        val wz = overlap1D(power.n, pg.power.n,
                            prod.origin.k - origin.k,
                            prod.origin.k - pg.origin.k, prod.exponent);

        return (Math.pow(Math.PI / prod.exponent, 1.5)
                             * Math.exp((-exponent * pg.exponent * radiusABSquared) / prod.exponent)
                             * wx * wy * wz);
    }


    public def mul(pg:PrimitiveGaussian) : PrimitiveGaussian {
        val gamma = exponent + pg.exponent;
        val newOrigin = Point3d(
                         (exponent * origin.i + pg.exponent * pg.origin.i) / gamma,
                         (exponent * origin.j + pg.exponent * pg.origin.j) / gamma,
                         (exponent * origin.k + pg.exponent * pg.origin.k) / gamma
                        );
        val pgres = new PrimitiveGaussian(newOrigin, Power(0,0,0), gamma, 0.0, false);

        return pgres;
    }

    /**
     * 1D overlap.
     *
     * <i> Taken from THO eq. 2.12 <i>
     */
    public def overlap1D(l1:Int, l2:Int, pax:Double, pbx:Double, gamma:Double) : Double {
        var sum:Double = 0.0;
        val k:Double = 1 + (Math.floor(0.5 * (l1 + l2)) as Int);

        // TODO: x10 parallel
        for(var i:Int = 0; i < k; i++) 
            sum += (MathUtil.binomialPrefactor(2 * i, l1, l2, pax, pbx)
                    * MathUtil.factorial2(2*i - 1)) / Math.pow(2 * gamma, i);

        return sum;
    }

    val PI_RAISE_TO_1DOT5 = Math.pow(Math.PI, 1.5);

    /**
     * The Kinetic Energy (KE) componant
     *
     * <i> Taken from THO eq. 2.14 <i>
     */
    public def kinetic(pg:PrimitiveGaussian) : Double {
        val l1 = power.l;
        val m1 = power.m;
        val n1 = power.n;
        val l2 = pg.power.l;
        val m2 = pg.power.m;
        val n2 = pg.power.n;

        var term:Double = pg.exponent * (2 * (l2+m2+n2) + 3) * ovrlp(pg);

        val origin = pg.origin;

        val p1 = new PrimitiveGaussian(origin, Power(l2+2, m2, n2), pg.exponent, pg.coefficient, false);
        val p2 = new PrimitiveGaussian(origin, Power(l2, m2+2, n2), pg.exponent, pg.coefficient, false);
        val p3 = new PrimitiveGaussian(origin, Power(l2, m2, n2+2), pg.exponent, pg.coefficient, false);

        term += -2.0 * Math.pow(pg.exponent, 2.0)
                     * (ovrlp(p1) + ovrlp(p2) + ovrlp(p3));

        val p4 = new PrimitiveGaussian(origin, Power(l2-2, m2, n2), pg.exponent, pg.coefficient, false);
        val p5 = new PrimitiveGaussian(origin, Power(l2, m2-2, n2), pg.exponent, pg.coefficient, false);
        val p6 = new PrimitiveGaussian(origin, Power(l2, m2, n2-2), pg.exponent, pg.coefficient, false);

        term += -0.5 * ((l2 * (l2 - 1)) * ovrlp(p4)
                        + (m2 * (m2 - 1)) * ovrlp(p5)
                        + (n2 * (n2 - 1)) * ovrlp(p6));

        return normalization * pg.normalization * term;
    }

    /**
     * The nuclear attraction term.
     *
     * <i> Taken from THO eq. 2.16 <i>
     */
    public def nuclear(pg:PrimitiveGaussian, centre:Point3d) :Double {
        val prod = mul(pg); 
        val rABSquared = origin.distanceSquared(pg.origin);
        val rCPSquared = centre.distanceSquared(prod.origin);

        val ax = new Array[Double](power.l + pg.power.l + 1);
        fillAArray(ax, power.l, pg.power.l,
                     prod.origin.i - origin.i,
                     prod.origin.i - pg.origin.i,
                     prod.origin.i - centre.i, prod.exponent);

        val ay = new Array[Double](power.m + pg.power.m + 1);
        fillAArray(ay, power.m, pg.power.m,
                     prod.origin.j - origin.j,
                     prod.origin.j - pg.origin.j,
                     prod.origin.j - centre.j, prod.exponent);

        val az = new Array[Double](power.n + pg.power.n + 1);
        fillAArray(az, power.n, pg.power.n,
                     prod.origin.k - origin.k,
                     prod.origin.k - pg.origin.k,
                     prod.origin.k - centre.k, prod.exponent);

        var sum:Double = 0.0;

        // TODO : x10 - parallel
        for(var i : Int = 0; i<ax.size; i++) {
            for(var j : Int = 0; j<ay.size; j++) {
                for(var k : Int = 0; k<az.size; k++) {
                    sum += ax(i) * ay(j) * az(k)
                          * IntegralsUtils.computeFGamma(i + j + k,
                                                        rCPSquared * prod.exponent);
                } // end for
            } // end for
        } // end for

        return (-normalization * pg.normalization * 2.0 * Math.PI / prod.exponent
                 * Math.exp(-exponent * pg.exponent * rABSquared / prod.exponent) * sum);
    }

    /**
     * <i> "THO eq. 2.18 and 3.1 <i>
     * note: assumes array a is already filled with zeros
     */
    private def fillAArray(a:Rail[Double], 
                            l1:Int, l2:Int, pa:Double, pb:Double,
                            cp:Double, gamma:Double) {

        // TODO : x10 - parallel
        for(var i:Int=0; i<a.size; i++) {
            for(var r:Int=0; r<((Math.floor(i/2.0)+1.0) as Int); r++) {
                for(var u:Int=0; u<((Math.floor((i-2.0 * r) / 2.0)+1.0) as Int); u++) {
                    a(i-2 * r-u) += Math.pow(-1.0, i) 
                         * MathUtil.binomialPrefactor(i, l1, l2, pa, pb)
                         * Math.pow(-1.0, u) * MathUtil.factorial(i)
                         * Math.pow(cp, i-2 * r-2 * u)
                         * Math.pow(0.25 / gamma, r + u)
                         / MathUtil.factorial(r)
                         / MathUtil.factorial(u)
                         / MathUtil.factorial(i-2 * r-2 * u);
                }
            }
        }
    }
}

