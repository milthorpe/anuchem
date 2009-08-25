/**
 * TwoElectronIntegrals.x10 
 * 
 * Evaluate 2E integrals 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class TwoElectronIntegrals { 
    var basisFunctions:BasisFunctions;
    var twoEInts:Array[Double]{rank==1};
    var direct:Boolean;
    var noOfIntegrals:Int;
    var contractedList:ArrayList[ContractedGaussian];

    public def this() {
        // Dummy constructor
        twoEInts = Array.make[Double]([0..1]);
        contractedList = new ArrayList[ContractedGaussian]();
        basisFunctions = null;
        direct = false;
        noOfIntegrals = 0;
    }

    public def make(bfs:BasisFunctions) { 
        basisFunctions = bfs;
        direct = false;

        val noOfBasisFunctions = basisFunctions.getBasisFunctions().size();
        noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;

        this.contractedList = null;

        // allocate required memory
        twoEInts = Array.make[Double]([0..noOfIntegrals]);
        compute2E();
    }

    public def make(bfs:BasisFunctions, isDirect:Boolean) {
        basisFunctions = bfs;
        direct = isDirect;

        val noOfBasisFunctions = basisFunctions.getBasisFunctions().size();
        noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;

        if (!direct) { 
           // allocate required memory
           twoEInts = Array.make[Double]([0..noOfIntegrals]);
           this.contractedList = null;
           compute2E();
        } else { 
           this.contractedList = basisFunctions.getBasisFunctions();
           twoEInts = Array.make[Double]([0..1]);
        } // end if
    }

    public def isDirect() : Boolean = direct;
    public def getNumberOfIntegrals() : Int = noOfIntegrals;

    public def compute2E(i:Int, j:Int, k:Int, l:Int) : Double {
        return coulomb(contractedList.get(i), contractedList.get(j), 
                       contractedList.get(k), contractedList.get(l));
    }

    protected def compute2E() : void {
        val bfs = basisFunctions.getBasisFunctions();
        val noOfBasisFunctions = bfs.size();
        

        var i:Int, j:Int, k:Int, l:Int, ij:Int, kl:Int, ijkl:Int;        
        
        var bfi:ContractedGaussian, bfj:ContractedGaussian, 
            bfk:ContractedGaussian, bfl:ContractedGaussian;
        
        // TODO: x10 - parallel
        // we only need i <= j, k <= l, and ij >= kl
        for(i=0; i<noOfBasisFunctions; i++) {
            bfi = bfs.get(i);
            
            for(j=0; j<(i+1); j++) {
                bfj = bfs.get(j);
                ij = i * (i+1) / 2+j;
                
                for(k=0; k<noOfBasisFunctions; k++) {
                    bfk = bfs.get(k);
                    
                    for(l=0; l<(k+1); l++) {
                        bfl = bfs.get(l);
                        
                        kl = k * (k+1) / 2+l;
                        if (ij >= kl) {
                            ijkl = IntegralsUtils.ijkl2intindex(i, j, k, l);
                            
                            // record the 2E integrals .. TODO
                            twoEInts(ijkl) = coulomb(bfi, bfj, bfk, bfl);
                        } // end if
                    } // end l loop
                } // end k loop
            } // end of j loop
        } // end of i loop       
    }

    public def getTwoElectronIntegrals() : Array[Double]{rank==1} = twoEInts;


    private def coulomb(a:ContractedGaussian, b:ContractedGaussian,
                        c:ContractedGaussian, d:ContractedGaussian) : Double {
         var jij:Double = 0.0;

         val aExps   = a.getExponents();
         val aCoefs  = a.getCoefficients();
         val aNorms  = a.getPrimNorms();
         val aOrigin = a.getOrigin();
         val aPower  = a.getPower();

         val bExps   = b.getExponents();
         val bCoefs  = b.getCoefficients();
         val bNorms  = b.getPrimNorms();
         val bOrigin = b.getOrigin();
         val bPower  = b.getPower();

         val cExps   = c.getExponents();
         val cCoefs  = c.getCoefficients();
         val cNorms  = c.getPrimNorms();
         val cOrigin = c.getOrigin();
         val cPower  = c.getPower();

         val dExps   = d.getExponents();
         val dCoefs  = d.getCoefficients();
         val dNorms  = d.getPrimNorms();
         val dOrigin = d.getOrigin();
         val dPower  = d.getPower();

         var i:Int, j:Int, k:Int, l:Int;
         var iaExp:Double, iaCoef:Double, iaNorm:Double,
             jbExp:Double, jbCoef:Double, jbNorm:Double,
             kcExp:Double, kcCoef:Double, kcNorm:Double;

         val na = aExps.size();
         val nb = bExps.size();
         val nc = cExps.size();
         val nd = dExps.size();

         // TODO: x10 parallel
         for(i=0; i<na; i++) {
             iaCoef = aCoefs.get(i);
             iaExp  = aExps.get(i);
             iaNorm = aNorms.get(i);

             for(j=0; j<nb; j++) {
                jbCoef = bCoefs.get(j);
                jbExp  = bExps.get(j);
                jbNorm = bNorms.get(j);

                for(k=0; k<nc; k++) {
                    kcCoef = cCoefs.get(k);
                    kcExp  = cExps.get(k);
                    kcNorm = cNorms.get(k);

                    for(l=0; l<nd; l++) {
                        jij += iaCoef * jbCoef * kcCoef
                               * dCoefs.get(l)
                               * coulombRepulsion(
                                         aOrigin, iaNorm, aPower, iaExp,
                                         bOrigin, jbNorm, bPower, jbExp,
                                         cOrigin, kcNorm, cPower, kcExp,
                                         dOrigin,
                                         dNorms.get(l),
                                         dPower,
                                         dExps.get(l)
                                        ); 
                    } // end l loop
                } // end k loop
             } // end j loop
         } // end i loop

         return (a.getNormalization() * b.getNormalization()
                 * c.getNormalization() * d.getNormalization() * jij);
    }

    private def coulombRepulsion(
                    a:Atom, aNorm:Double, aPower:Power, aAlpha:Double,
                    b:Atom, bNorm:Double, bPower:Power, bAlpha:Double,
                    c:Atom, cNorm:Double, cPower:Power, cAlpha:Double,
                    d:Atom, dNorm:Double, dPower:Power, dAlpha:Double) : Double {

        val radiusABSquared = a.distanceSquaredFrom(b);
        val radiusCDSquared = c.distanceSquaredFrom(d);

        val p = IntegralsUtils.gaussianProductCenter(aAlpha, a, bAlpha, b);
        val q = IntegralsUtils.gaussianProductCenter(cAlpha, c, dAlpha, d);

        val radiusPQSquared = p.distanceSquaredFrom(q);

        val gamma1 = aAlpha + bAlpha;
        val gamma2 = cAlpha + dAlpha;
        val delta  = 0.25 * (1/gamma1 + 1/gamma2);

        val bx = constructBArray(
                   aPower.getL(), bPower.getL(), cPower.getL(), dPower.getL(),
                   p.getX(), a.getX(), b.getX(), q.getX(), c.getX(), d.getX(),
                   gamma1, gamma2, delta);

        val by = constructBArray(
                   aPower.getM(), bPower.getM(), cPower.getM(), dPower.getM(),
                   p.getY(), a.getY(), b.getY(), q.getY(), c.getY(), d.getY(),
                   gamma1, gamma2, delta);

        val bz = constructBArray(
                   aPower.getN(), bPower.getN(), cPower.getN(), dPower.getN(),
                   p.getZ(), a.getZ(), b.getZ(), q.getZ(), c.getZ(), d.getZ(),
                   gamma1, gamma2, delta);

        var sum:Double = 0.0;
        var i:Int, j:Int, k:Int;
        val nbx = bx.region.max(0);
        val nby = by.region.max(0);
        val nbz = bz.region.max(0);
        // TODO: x10 parallel
        for(i=0; i<nbx; i++) {
            for(j=0; j<nby; j++) {
                for(k=0; k<nbz; k++) {
                    sum += bx(i) * by(j) * bz(k)
                           * IntegralsUtils.computeFGamma(
                                            i+j+k, 0.25*radiusPQSquared/delta);
                } // end for
            } // end for
        } // end for

        return (2 * Math.pow(Math.PI, 2.5)
                  / (gamma1 * gamma2 * Math.sqrt(gamma1+gamma2))
                  * Math.exp(-aAlpha*bAlpha*radiusABSquared/gamma1)
                  * Math.exp(-cAlpha*dAlpha*radiusCDSquared/gamma2)
                  * sum * aNorm * bNorm * cNorm * dNorm);
    }

    /**
     * Construct B array.
     *
     * <i> THO eq. 2.22 </i>
     */
    private def constructBArray(l1:Int, l2:Int, l3:Int, l4:Int,
                  p:Double, a:Double, b:Double, q:Double, c:Double, d:Double,
                  g1:Double, g2:Double, delta:Double) : Array[Double]{rank==1} {
        val iMax = l1+l2+l3+l4+1;
        val bArr:Array[Double]{rank==1} = Array.make[Double]([0..iMax]);

        var i1:Int, i2:Int, r1:Int, r2:Int, u:Int, index:Int;

        // TODO: x10 parallel
        for(i1=0; i1<(l1+l2+1); i1++) {
            for(i2=0; i2<(l3+l4+1); i2++) {
                for(r1=0; r1<(i1/2+1); r1++) {
                    for(r2=0; r2<(i2/2+1); r2++) {
                        for(u=0; u<((i1+i2)/2-r1-r2+1); u++) {
                            index = i1 + i2 - 2*(r1+r2) - u;

                            bArr(index) += constructBTerm(i1, i2, r1, r2, u,
                                              l1, l2, l3, l4, p, a, b, q, c, d,
                                              g1, g2, delta);
                        } // end u loop
                    } // end r2 loop
                } // end r1 loop
            } // end i2 loop
        } // end i1 loop

        return bArr;
    }

    /**
     * Construct the B term
     *
     * <i> THO eq. 2.22 </i>
     */
    private def constructBTerm(i1:Int, i2:Int, r1:Int, r2:Int, u:Int,
         l1:Int, l2:Int, l3:Int, l4:Int, px:Double, ax:Double, bx:Double,
         qx:Double, cx:Double, dx:Double, gamma1:Double, gamma2:Double,
         delta:Double) : Double {

        return (functionB(i1, l1, l2, px, ax, bx, r1, gamma1)
                * Math.pow(-1, i2)
                * functionB(i2, l3, l4, qx, cx, dx, r2, gamma2)
                * Math.pow(-1, u)
                * MathUtil.factorialRatioSquared(i1+i2-2*(r1+r2), u)
                * Math.pow(qx-px, i1+i2-2*(r1+r2)-2*u)
                / Math.pow(delta, i1+i2-2*(r1+r2)-u));
    }

    /**
     * the function B
     */
    private def functionB(i:Int, l1:Int, l2:Int, p:Double, a:Double,
                          b:Double, r:Int, g:Double) : Double {
        return (MathUtil.binomialPrefactor(i, l1, l2, p-a, p-b)
                * functionB0(i, r, g));
    }

    /**
     * the function B0
     */
    private def functionB0(i:Int, r:Int, g:Double) : Double {
        return (MathUtil.factorialRatioSquared(i, r)
                * Math.pow(4*g, r-i));
    }
}

