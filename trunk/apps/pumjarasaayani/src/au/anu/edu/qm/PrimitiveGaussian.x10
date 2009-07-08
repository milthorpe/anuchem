package au.anu.edu.qm;

public class PrimitiveGaussian { 
    var origin:Atom;
    var power:Power;
    var exponent:Double;
    var coefficient:Double;
    var normalization:Double;

    public def this(atm:Atom, pwr:Power, exp:Double, coeff:Double) { 
        origin = atm;
        power  = pwr;
        exponent    = exp;
        coefficient = coeff;

        normalize();
    } 

    public def getOrigin() : Atom = origin;
    public def getPower()  : Power = power;
    public def getExponent() : Double = exponent;
    public def getCoefficient() : Double = coefficient;
    public def getNormalization() : Double = normalization;
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();


    public def overlap(pg:PrimitiveGaussian) : Double {
        return normalization * pg.normalization * ovrlp(pg);
    }

    private def ovrlp(pg:PrimitiveGaussian) : Double {
        val radiusABSquared = origin.distanceSquaredFrom(pg.origin);
        val prod = mul(pg);

        val wx = overlap1D(power.getL(), pg.power.getL(),
                            prod.origin.getX() - origin.getX(),
                            prod.origin.getX() - pg.origin.getX(), prod.exponent);
        val wy = overlap1D(power.getM(), pg.power.getM(),
                            prod.origin.getY() - origin.getY(),
                            prod.origin.getY() - pg.origin.getY(), prod.exponent);
        val wz = overlap1D(power.getN(), pg.power.getN(),
                            prod.origin.getZ() - origin.getZ(),
                            prod.origin.getZ() - pg.origin.getZ(), prod.exponent);

        return (Math.pow(Math.PI / prod.exponent, 1.5)
                             * Math.exp((-exponent * pg.exponent * radiusABSquared) / prod.exponent)
                             * wx * wy * wz);
    }


    public def mul(pg:PrimitiveGaussian) : PrimitiveGaussian {
        val gamma = exponent + pg.exponent;
        val newOrigin = new Atom(
                         (exponent * origin.getX() + pg.exponent * pg.origin.getX()) / gamma,
                         (exponent * origin.getY() + pg.exponent * pg.origin.getY()) / gamma,
                         (exponent * origin.getZ() + pg.exponent * pg.origin.getZ()) / gamma
                        );
        return new PrimitiveGaussian(newOrigin, new Power(0,0,0), gamma, 0.0);
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
    public def kinetic(pg:PrimitiveGaussian) : Double {
        val l1 = power.getL();
        val m1 = power.getM();
        val n1 = power.getN();
        val l2 = pg.power.getL();
        val m2 = pg.power.getM();
        val n2 = pg.power.getN();

        var term:Double = pg.exponent * (2 * (l2+m2+n2) + 3) * ovrlp(pg);

        term += -2.0 * Math.pow(pg.exponent, 2.0)
                     * (ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2+2, m2, n2), 
                                                      pg.exponent, pg.coefficient))
                        + ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2, m2+2, n2),   
                                                      pg.exponent, pg.coefficient))
                        + ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2, m2, n2+2),   
                                                      pg.exponent, pg.coefficient)));
        term += -0.5 * ((l2 * (l2 - 1)) * ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2-2, m2, n2),   
                                                  pg.exponent, pg.coefficient))
                        + (m2 * (m2 - 1)) * ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2, m2-2, n2),   
                                                      pg.exponent, pg.coefficient)) 
                        + (n2 * (n2 - 1)) * ovrlp(new PrimitiveGaussian(pg.origin, new Power(l2, m2, n2-2),   
                                                      pg.exponent, pg.coefficient)));

        return normalization * pg.normalization * term;
    }

    /**
     * The nuclear attraction term.
     *
     * <i> Taken from THO eq. 2.12 <i>
     */
    public def nuclear(pg:PrimitiveGaussian, center:Atom) :Double {
        val prod = mul(pg); 
        val rABSquared = origin.distanceSquaredFrom(pg.origin);
        val rCPSquared = center.distanceSquaredFrom(prod.origin);
        var nterm:Double = 0.0;

        val ax = constructAArray(power.getL(), pg.power.getL(),
                                 prod.origin.getX() - origin.getX(),
                                 prod.origin.getX() - pg.origin.getX(),
                                 prod.origin.getX() - center.getX(), prod.exponent);

        val ay = constructAArray(power.getM(), pg.power.getM(),
                                 prod.origin.getY() - origin.getY(),
                                 prod.origin.getY() - pg.origin.getY(),
                                 prod.origin.getY() - center.getY(), prod.exponent);

        val az = constructAArray(power.getN(), pg.power.getN(),
                                 prod.origin.getZ() - origin.getZ(),
                                 prod.origin.getZ() - pg.origin.getZ(),
                                 prod.origin.getZ() - center.getZ(), prod.exponent);

        var sum:Double = 0.0;
        var i:Int, j:Int, k:Int;
        val nax = ax.region.max(0);
        val nay = ay.region.max(0);
        val naz = az.region.max(0);

        // TODO : x10 - parallel
        for(i = 0; i<nax; i++) {
            for(j = 0; j<nay; j++) {
                for(k = 0; k<naz; k++) {
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
                                cp:Double, gamma:Double) : Array[Double]{rank==1} {
        val iMax    = l1 + l2 + 1;
        val a:Array[Double]{rank==1} = Array.make[Double]([0..iMax]);

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

