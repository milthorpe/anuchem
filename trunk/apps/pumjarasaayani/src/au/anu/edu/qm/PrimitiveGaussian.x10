/**
 * PrimitiveGaussian.x10
 *
 * Represents a primitiv gaussian function 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.vector.Point3d;

public class PrimitiveGaussian { 
    global val origin:Point3d;
    global val power:Power;
    global val exponent:Double;
    global var coefficient:Double;
    global var normalization:Double;

    public def this(origin:Point3d, power:Power, exponent:Double, coefficient:Double) { 
        this.origin = origin;
        this.power = power;
        this.exponent = exponent;
        this.coefficient = coefficient;
    } 

    public def getOrigin() = origin;
    public def getPower() = power;
    public def getExponent() = exponent;
    public def getCoefficient() = coefficient;
    public def getNormalization() = normalization;
    public def setNormalization(n:Double) : Void { normalization = n; }
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    public def setCoefficient(coef:Double) { coefficient = coef; }

    public def overlap(pg:PrimitiveGaussian{self.at(this)}) : Double {
        return normalization * pg.normalization * ovrlp(pg);
    }

    private def ovrlp(pg:PrimitiveGaussian{self.at(this)}) : Double {
        val radiusABSquared = origin.distanceSquared(pg.origin);
        val prod = mul(pg);

        val wx = overlap1D(power.getL(), pg.power.getL(),
                            prod.origin.i - origin.i,
                            prod.origin.i - pg.origin.i, prod.exponent);
        val wy = overlap1D(power.getM(), pg.power.getM(),
                            prod.origin.j - origin.j,
                            prod.origin.j - pg.origin.j, prod.exponent);
        val wz = overlap1D(power.getN(), pg.power.getN(),
                            prod.origin.k - origin.k,
                            prod.origin.k - pg.origin.k, prod.exponent);

        return (Math.pow(Math.PI / prod.exponent, 1.5)
                             * Math.exp((-exponent * pg.exponent * radiusABSquared) / prod.exponent)
                             * wx * wy * wz);
    }


    public def mul(pg:PrimitiveGaussian{self.at(this)}) : PrimitiveGaussian{self.at(this)} {
        val gamma = exponent + pg.exponent;
        val newOrigin = new Point3d(
                         (exponent * origin.i + pg.exponent * pg.origin.i) / gamma,
                         (exponent * origin.j + pg.exponent * pg.origin.j) / gamma,
                         (exponent * origin.k + pg.exponent * pg.origin.k) / gamma
                        );
        val pgres:PrimitiveGaussian{self.at(this)} = new PrimitiveGaussian(newOrigin, new Power(0,0,0), gamma, 0.0);

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

    public def normalize() : void {
       val l = power.getL(); 
       val m = power.getM();
       val n = power.getN();       
        
       normalization = Math.sqrt(Math.pow(2, 2 * (l + m + n) + 1.5) *
                                   Math.pow(exponent, l + m + n + 1.5) /
                                   MathUtil.factorial2(2 * l - 1) / 
                                   MathUtil.factorial2(2 * m - 1) /
                                   MathUtil.factorial2(2 * n - 1) /
                                   PI_RAISE_TO_1DOT5);
    }

    /**
     * The Kinetic Energy (KE) componant
     *
     * <i> Taken from THO eq. 2.12 <i>
     */
    public def kinetic(pg:PrimitiveGaussian{self.at(this)}) : Double {
        val l1 = power.getL();
        val m1 = power.getM();
        val n1 = power.getN();
        val l2 = pg.power.getL();
        val m2 = pg.power.getM();
        val n2 = pg.power.getN();

        var term:Double = pg.exponent * (2 * (l2+m2+n2) + 3) * ovrlp(pg);

        val origin = pg.origin;

        val p1 = new PrimitiveGaussian(origin, new Power(l2+2, m2, n2), pg.exponent, pg.coefficient);
        val p2 = new PrimitiveGaussian(origin, new Power(l2, m2+2, n2), pg.exponent, pg.coefficient);
        val p3 = new PrimitiveGaussian(origin, new Power(l2, m2, n2+2), pg.exponent, pg.coefficient);

        term += -2.0 * Math.pow(pg.exponent, 2.0)
                     * (ovrlp(p1) + ovrlp(p2) + ovrlp(p3));

        val p4 = new PrimitiveGaussian(origin, new Power(l2-2, m2, n2), pg.exponent, pg.coefficient);
        val p5 = new PrimitiveGaussian(origin, new Power(l2, m2-2, n2), pg.exponent, pg.coefficient);
        val p6 = new PrimitiveGaussian(origin, new Power(l2, m2, n2-2), pg.exponent, pg.coefficient);

        term += -0.5 * ((l2 * (l2 - 1)) * ovrlp(p4)
                        + (m2 * (m2 - 1)) * ovrlp(p5)
                        + (n2 * (n2 - 1)) * ovrlp(p6));

        return normalization * pg.normalization * term;
    }

    /**
     * The nuclear attraction term.
     *
     * <i> Taken from THO eq. 2.12 <i>
     */
    public def nuclear(pg:PrimitiveGaussian{self.at(this)}, center:Point3d) :Double {
        val prod = mul(pg); 
        val rABSquared = origin.distanceSquared(pg.origin);
        val rCPSquared = center.distanceSquared(prod.origin);
        var nterm:Double = 0.0;

        val ax = constructAArray(power.getL(), pg.power.getL(),
                                 prod.origin.i - origin.i,
                                 prod.origin.i - pg.origin.i,
                                 prod.origin.i - center.i, prod.exponent);

        val ay = constructAArray(power.getM(), pg.power.getM(),
                                 prod.origin.j - origin.j,
                                 prod.origin.j - pg.origin.j,
                                 prod.origin.j - center.j, prod.exponent);

        val az = constructAArray(power.getN(), pg.power.getN(),
                                 prod.origin.k - origin.k,
                                 prod.origin.k - pg.origin.k,
                                 prod.origin.k - center.k, prod.exponent);

        var sum:Double = 0.0;

        // TODO : x10 - parallel
        for(var i : Int = 0; i<ax.length; i++) {
            for(var j : Int = 0; j<ay.length; j++) {
                for(var k : Int = 0; k<az.length; k++) {
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
     */
    private def constructAArray(l1:Int, l2:Int, pa:Double, pb:Double,
                                cp:Double, gamma:Double) {
        val iMax    = l1 + l2 + 1;
        val a = Rail.make[Double](iMax);

        var i:Int, r:Int, u:Int, index:Int;

        for(i=0; i<iMax; i++) a(i) = 0.0;

        // TODO : x10 - parallel
        for(i=0; i<iMax; i++) {
            for(r=0; r<((Math.floor(i/2.0)+1.0) as Int); r++) {
                for(u=0; u<((Math.floor((i-2.0 * r) / 2.0)+1.0) as Int); u++) {
                    index = i-2 * r-u;

                    a(index) += constructATerm(i, r, u, l1, l2,
                                               pa, pb, cp, gamma);
                } // end for
            } // end for
        } // end for

        return a;
    }

    /**
     * the A term <br>
     * <i> "THO eq. 2.18 <i>
     */
    private def constructATerm(i:Int, r:Int, u:Int, l1:Int, l2:Int,
                               pax:Double, pbx:Double, cpx:Double,
                               gamma:Double) : Double {
        return (Math.pow(-1.0, i) * MathUtil.binomialPrefactor(i, l1, l2, pax, pbx)
                 * Math.pow(-1.0, u) * MathUtil.factorial(i)
                 * Math.pow(cpx, i-2 * r-2 * u)
                 * Math.pow(0.25 / gamma, r + u)
                 / MathUtil.factorial(r)
                 / MathUtil.factorial(u)
                 / MathUtil.factorial(i-2 * r-2 * u));
    }
}

