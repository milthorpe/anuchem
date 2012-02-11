/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 * (C) Copyright Josh Milthorpe 2011-2012.
 */
package au.edu.anu.qm;

import x10.compiler.Inline;
import x10x.vector.Point3d;

/**
 * Represents a primitive gaussian function
 *
 * @author: V.Ganesh, milthorpe
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
            this.normalization = PrimitiveGaussian.getNormalization(power, exponent);
        } else {
            this.normalization = 0.0; // TODO - should be 1.0?
        }
    }


    static @Inline def getNormalization(power:Power, exponent:Double):Double {
        return Math.sqrt(Math.pow(2, 2 * (power.l + power.m + power.n) + 1.5) *
                           Math.pow(exponent, power.l + power.m + power.n + 1.5) /
                           MathUtil.factorial2(2 * power.l - 1) / 
                           MathUtil.factorial2(2 * power.m - 1) /
                           MathUtil.factorial2(2 * power.n - 1) /
                           Math.pow(Math.PI, 1.5));
    }

    public def overlap(pg:PrimitiveGaussian) : Double {
        return normalization * pg.normalization * PrimitiveGaussian.overlap(exponent, origin, power, pg.exponent, pg.origin, pg.power);
    }

    static @Inline def overlap(expI:Double, originI:Point3d, powerI:Power, expJ:Double, originJ:Point3d, powerJ:Power) {
        val radiusABSquared = originI.distanceSquared(originJ);
        val gamma = expI + expJ;
        val newOrigin = Point3d(
                 (expI * originI.i + expJ * originJ.i) / gamma,
                 (expI * originI.j + expJ * originJ.j) / gamma,
                 (expI * originI.k + expJ * originJ.k) / gamma
        );

        val wx = PrimitiveGaussian.overlap1D(powerI.l, powerJ.l,
                            newOrigin.i - originI.i,
                            newOrigin.i - originJ.i, gamma);
        val wy = PrimitiveGaussian.overlap1D(powerI.m, powerJ.m,
                            newOrigin.j - originI.j,
                            newOrigin.j - originJ.j, gamma);
        val wz = PrimitiveGaussian.overlap1D(powerI.n, powerJ.n,
                            newOrigin.k - originI.k,
                            newOrigin.k - originJ.k, gamma);

        return Math.pow(Math.PI / gamma, 1.5)
                     * Math.exp((-expI * expJ * radiusABSquared) / gamma)
                     * wx * wy * wz;
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
    public static def overlap1D(l1:Int, l2:Int, pax:Double, pbx:Double, gamma:Double) : Double {
        var sum:Double = 0.0;
        val k:Double = 1.0 + Math.floor(0.5 * (l1 + l2));

        // TODO: x10 parallel
        for(var i:Int = 0; i < k; i++) 
            sum += (MathUtil.binomialPrefactor(2 * i, l1, l2, pax, pbx)
                    * MathUtil.factorial2(2*i - 1)) / Math.pow(2 * gamma, i);

        return sum;
    }

    public def kinetic(pg:PrimitiveGaussian) : Double {
        return PrimitiveGaussian.kinetic(exponent, origin, power, normalization, pg.exponent, pg.origin, pg.power, pg.normalization);
    }

    /**
     * The Kinetic Energy (KE) componant
     *
     * <i> Taken from THO eq. 2.14 <i>
     */
    static @Inline def kinetic(expI:Double, originI:Point3d, powerI:Power, normI:Double, expJ:Double, originJ:Point3d, powerJ:Power, normJ:Double) {
        val l2 = powerJ.l;
        val m2 = powerJ.m;
        val n2 = powerJ.n;

        var term:Double = expJ * (2 * (l2+m2+n2) + 3) * overlap(expI, originI, powerI, expJ, originJ, powerJ);
        term += -2.0 * expJ*expJ
                     * (overlap(expI, originI, powerI, expJ, originJ, Power(l2+2, m2, n2))
                      + overlap(expI, originI, powerI, expJ, originJ, Power(l2, m2+2, n2))
                      + overlap(expI, originI, powerI, expJ, originJ, Power(l2, m2, n2+2)));

        term += -0.5 * ((l2 * (l2-1)) * overlap(expI, originI, powerI, expJ, originJ, Power(l2-2, m2, n2))
                      + (m2 * (m2-1)) * overlap(expI, originI, powerI, expJ, originJ, Power(l2, m2-2, n2))
                      + (n2 * (n2-1)) * overlap(expI, originI, powerI, expJ, originJ, Power(l2, m2, n2-2)));

        return normI * normJ * term;
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
            for(var r:Int=0; r<(i/2+1); r++) {
                for(var u:Int=0; u<((i-2*r)/2+1); u++) {
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

