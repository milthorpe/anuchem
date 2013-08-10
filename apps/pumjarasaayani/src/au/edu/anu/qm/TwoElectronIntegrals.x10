/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.compiler.Inline;
import x10.regionarray.Array;
import x10.regionarray.Region;
import x10.util.ArrayList;

import x10.matrix.DenseMatrix;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.Molecule;

/**
 * This class calculates two-electron integrals
 *
 * @author: V.Ganesh
 *
 * Old code based on original integral evaluation paper:
 *   H. Taketa, S. Huzinaga, and K. O-ohata. H. (1966).
 *   "Gaussian-Expansion Methods for Molecular Integrals".
 *   Phys. Soc. Japan, 21 (11) p. 2313.
 *
 * New Murchie-Davidson (MD) code based on:
 *   McMurchie, L. E. and Davidson, E. R. (1978).
 *   "One- and two-electron integrals over cartesian gaussian functions".
 *   J. Comp. Phys, 26 (2) pp. 218-231.
 */
public class TwoElectronIntegrals {
    private static SQ2PI = 0.797884560802865; //Math.pow((2.0/Math.PI), 0.5);

    private val gmt:Rail[Double], zeroM:Rail[Double];

    private val rM:Array[Double](2){rect,zeroBased};
    private val pqInts:Array[Double](4){rect,zeroBased};
    private val npint:Array[Double](2){rect,zeroBased};
    private val pcdint:Array[Double](4){rect,zeroBased};

    private val maxam:Long, maxam2:Long, maxam2M:Long, maxam2N:Long;

    private val normFactors:Rail[Double];
    private val THRESH:Double;
    private val TCrit:Double;

    private val omega:Double;

    /**
     * @param maxam maximum angular momentum (determines total number of integrals)
     * @param normFactors normalization factors for integrals of different angular momenta
     */
    public def this(maxam:Int, normFactors:Rail[Double], omega:Double, Th:Double) {
        // allocate scratch memory
        this.omega=omega;
        this.maxam = maxam;
        maxam2 = 2n*maxam;
        val maxam4 = 4n*maxam;
        val maxamN = ((maxam+1)*(maxam+2)/2);
        maxam2M  = ((maxam2+1)*(maxam2+2)/2);
        val maxam2N  = ((maxam2+1)*(maxam2M+1));
        this.maxam2N = maxam2N;

        this.normFactors = normFactors;
        this.THRESH = Th;
        this.TCrit = - Math.log(Th);

        // Console.OUT.println("alloc: " + maxam + " " + maxam2N);

        gmt    = new Rail[Double](maxam4+1);
        zeroM  = new Rail[Double](maxam4+1);

        rM     = new Array[Double](Region.make(0..(maxam4+1), 0..((maxam4+1)*(maxam4+2)/2)));
        pqInts = new Array[Double](Region.make(0..maxam2, 0..maxam2M, 0..maxam2, 0..maxam2M));
        npint  = new Array[Double](Region.make(0..maxam2, 0..maxam2M));
        pcdint = new Array[Double](Region.make(0..(maxamN+1), 0..(maxamN+1), 0..maxam2, 0..maxam2M));

        // Console.OUT.println("alloc2: " + pcdint.region.size());
    }

    /**
     * Compute all two electron integrals for given shells and store in J,K matrices.
     * Uses MD recurrence relations to evaluate higher angular momentum integrals. 
     */
    public def compute2EAndRecord(a:ContractedGaussian, b:ContractedGaussian, 
                                  c:ContractedGaussian, d:ContractedGaussian, 
                                  shellList:ShellList, 
                                  jMat:DenseMatrix, kMat:DenseMatrix,
                                  dMat:Density):Int {
        val aCen  = a.origin;
        val bCen  = b.origin;
        val cCen  = c.origin;
        val dCen  = d.origin;

        val dAng  = d.getMaximumAngularMomentum();
        val cAng  = c.getMaximumAngularMomentum();
        val bAng  = b.getMaximumAngularMomentum();
        val aAng  = a.getMaximumAngularMomentum();

        val shellA = shellList.getPowers(aAng);
        val shellB = shellList.getPowers(bAng);
        val shellC = shellList.getPowers(cAng);
        val shellD = shellList.getPowers(dAng);

        val dLim = ((dAng+1n)*(dAng+2n)/2n);
        val cLim = ((cAng+1n)*(cAng+2n)/2n);
        val bLim = ((bAng+1n)*(bAng+2n)/2n);
        val aLim = ((aAng+1n)*(aAng+2n)/2n);

        val nTot = aLim*bLim*cLim*dLim;

        val dStrt = d.intIndex;
        val cStrt = c.intIndex;
        val bStrt = b.intIndex;
        val aStrt = a.intIndex;

        val angMomAB = aAng + bAng;
        val angMomCD = cAng + dAng;
        val angMomABCD = angMomAB+angMomCD;

        val radiusABSquared = a.distanceSquaredFrom(b); 
        val radiusCDSquared = c.distanceSquaredFrom(d);

        val aExps = a.exponents;
        val aCoeffs = a.coefficients;
        val bExps = b.exponents;
        val bCoeffs = b.coefficients;
        val cExps = c.exponents;
        val cCoeffs = c.coefficients;
        val dExps = d.exponents;
        val dCoeffs = d.coefficients;

        for(ap in 0..(aExps.size-1)) {
            val aAlpha = aExps(ap);
            val aCoeff = aCoeffs(ap);

           for(bp in 0..(bExps.size-1)) {
             pcdint.clear();
             val bAlpha = bExps(bp);
             val gamma1 = (aAlpha + bAlpha);
             val sigmaP = 1.0 / gamma1;
             val bCoeff = bCoeffs(bp);

             // val p = gaussianProductCentre(aAlpha, aCen, bAlpha, bCen);

             val p = Point3d(
                         (aAlpha * aCen.i + bAlpha * bCen.i) * sigmaP,
                         (aAlpha * aCen.j + bAlpha * bCen.j) * sigmaP,
                         (aAlpha * aCen.k + bAlpha * bCen.k) * sigmaP
                       );

             val Gab = Math.exp(-aAlpha*bAlpha*radiusABSquared*sigmaP);
             val pgx = Math.PI * sigmaP;
             val Up = aCoeff*bCoeff*Gab*Math.sqrt(pgx*pgx*pgx);
             // Console.OUT.println("Coeff: " + aCoeff + " " + bCoeff);
             // Console.OUT.println("sigmaP, Gab, Up: " + sigmaP + " " + Gab + " " + Up);

             for(cp in 0..(cExps.size-1)) {
               val cAlpha = cExps(cp);
               val cCoeff = cCoeffs(cp);

               for(dp in 0..(dExps.size-1)) {
                 val dAlpha = dExps(dp);
                 val dCoeff = dCoeffs(dp);

                 val gamma2 = (cAlpha + dAlpha);
                 val sigmaQ = 1.0 / gamma2;
                 val eta = (gamma1*gamma2)/(gamma1+gamma2);

                 // val q = gaussianProductCentre(cAlpha, cCen, dAlpha, dCen);
                 val q = Point3d(
                           (cAlpha * cCen.i + dAlpha * dCen.i) * sigmaQ,
                           (cAlpha * cCen.j + dAlpha * dCen.j) * sigmaQ,
                           (cAlpha * cCen.k + dAlpha * dCen.k) * sigmaQ
                         );

                 val Gcd = Math.exp(-cAlpha*dAlpha*radiusCDSquared*sigmaQ);
                 val qgx = Math.PI * sigmaQ;
                 val Uq  = cCoeff*dCoeff*Gcd*Math.sqrt(qgx*qgx*qgx); 
                 // Console.OUT.println("Coeff: " + cCoeff + " " + dCoeff);
                 // Console.OUT.println("sigmaQ, Gcd, Uq: " + sigmaQ + " " + Gcd + " " + Uq);

                 val r = q.vector(p);
                 val radiusPQSquared = r.lengthSquared();
                 val T = radiusPQSquared * eta;

                 // Console.OUT.println("Computing [0]m");
                 // Console.OUT.println("T: " + T + "Up, Uq, Upq: " + Up + " " + Uq + " " + Upq);

                 // compute [0]m
                 val Upq = Up*Uq;
                 computeZeroM(angMomABCD, T, Upq, eta);

                 // Console.OUT.println("Computing [r]m");
                 // Console.OUT.println("abcd_ang: " + angMomABCD);

                 // form [r]m using MD (recursion)
                 computeRm(angMomABCD, shellList, r);
 
                 // Console.OUT.println("Computing [p|q] ");

                 // form [p|q] 
                 computePq(angMomAB, angMomCD, shellList);

                 // Console.OUT.println("Computing [p|cd] ");

                 // form [p|cd]
                 computePcd(angMomAB, sigmaQ, q, dLim, cLim,
                            shellD, shellC, dCen, cCen);
               } // d
              } // c

              // form [ab|cd], 
              // Console.OUT.println("Computing [ab|cd] ");

              computeAbcd(dLim, cLim, bLim, aLim,
                          dStrt, cStrt, bStrt, aStrt,
                          shellB, shellA,
                          aCen, bCen, p, sigmaP,
                          jMat, kMat, dMat); 
           }
        }
        return nTot;
    }

    /** A modification of above function */
    public def compute2EAndRecord2(a:ContractedGaussian, b:ContractedGaussian, 
                                  c:ContractedGaussian, d:ContractedGaussian, 
                                  shellList:ShellList, 
                                  jMat:DenseMatrix, kMat:DenseMatrix,
                                  dMat:Density,
                                  radiusABSquared:Double, 
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int):Int {
        val aCen  = a.origin;
        val bCen  = b.origin;
        val cCen  = c.origin;
        val dCen  = d.origin;

        val shellA = shellList.getPowers(aAng);
        val shellB = shellList.getPowers(bAng);
        val shellC = shellList.getPowers(cAng);
        val shellD = shellList.getPowers(dAng);

        val dLim = ((dAng+1n)*(dAng+2n)/2n);
        val cLim = ((cAng+1n)*(cAng+2n)/2n);
        val nTot = aLim*bLim*cLim*dLim;

        val angMomCD = cAng + dAng;
        val angMomABCD = angMomAB+angMomCD;

        val radiusCDSquared = c.distanceSquaredFrom(d);

        val aExps = a.exponents;
        val aCoeffs = a.coefficients;
        val bExps = b.exponents;
        val bCoeffs = b.coefficients;
        val cExps = c.exponents;
        val cCoeffs = c.coefficients;
        val dExps = d.exponents;
        val dCoeffs = d.coefficients;

        for(ap in 0..(aExps.size-1)) {
          val aAlpha = aExps(ap);
          val aCoeff = aCoeffs(ap);

          for(bp in 0..(bExps.size-1)) {
             pcdint.clear();
             val bAlpha = bExps(bp);
             val gamma1 = (aAlpha + bAlpha);
             val sigmaP = 1.0 / gamma1;
             val bCoeff = bCoeffs(bp);

             // val p = gaussianProductCentre(aAlpha, aCen, bAlpha, bCen);

             val p = Point3d(
                         (aAlpha * aCen.i + bAlpha * bCen.i) * sigmaP,
                         (aAlpha * aCen.j + bAlpha * bCen.j) * sigmaP,
                         (aAlpha * aCen.k + bAlpha * bCen.k) * sigmaP
                       );

             val Gab = Math.exp(-aAlpha*bAlpha*radiusABSquared*sigmaP);
             val pgx = Math.PI * sigmaP;
             val Up = aCoeff*bCoeff*Gab*Math.sqrt(pgx*pgx*pgx);
             // Console.OUT.println("Coeff: " + aCoeff + " " + bCoeff);
             // Console.OUT.println("sigmaP, Gab, Up: " + sigmaP + " " + Gab + " " + Up);

             for(cp in 0..(cExps.size-1)) {
               val cAlpha = cExps(cp);
               val cCoeff = cCoeffs(cp);

               for(dp in 0..(dExps.size-1)) {
                 val dAlpha = dExps(dp);
                 val dCoeff = dCoeffs(dp);

                 val gamma2 = (cAlpha + dAlpha);
                 val sigmaQ = 1.0 / gamma2;
                 val eta = (gamma1*gamma2)/(gamma1+gamma2);

                 // val q = gaussianProductCentre(cAlpha, cCen, dAlpha, dCen);
                 val q = Point3d(
                           (cAlpha * cCen.i + dAlpha * dCen.i) * sigmaQ,
                           (cAlpha * cCen.j + dAlpha * dCen.j) * sigmaQ,
                           (cAlpha * cCen.k + dAlpha * dCen.k) * sigmaQ
                         );

                 val Gcd = Math.exp(-cAlpha*dAlpha*radiusCDSquared*sigmaQ);
                 val qgx = Math.PI * sigmaQ;
                 val Uq  = cCoeff*dCoeff*Gcd*Math.sqrt(qgx*qgx*qgx); 
                 // Console.OUT.println("Coeff: " + cCoeff + " " + dCoeff);
                 // Console.OUT.println("sigmaQ, Gcd, Uq: " + sigmaQ + " " + Gcd + " " + Uq);

                 val r = q.vector(p);
                 val radiusPQSquared = r.lengthSquared();  
                 val T = radiusPQSquared * eta;

                 // Console.OUT.println("Computing [0]m");
                 // Console.OUT.println("T: " + T + "Up, Uq, Upq: " + Up + " " + Uq + " " + Upq);

                 // compute [0]m
                 val Upq = Up*Uq;
                 computeZeroM(angMomABCD, T, Upq, eta);

                 // Console.OUT.println("Computing [r]m");
                 // Console.OUT.println("abcd_ang: " + angMomABCD);

                 // form [r]m using MD (recursion)
                 computeRm(angMomABCD, shellList, r);
 
                 // Console.OUT.println("Computing [p|q] ");

                 // form [p|q] 
                 computePq(angMomAB, angMomCD, shellList);

                 // Console.OUT.println("Computing [p|cd] ");

                 // form [p|cd]
                 computePcd(angMomAB, sigmaQ, q, dLim, cLim,
                            shellD, shellC, dCen, cCen);
               } // d
              } // c

              // form [ab|cd], 
              // Console.OUT.println("Computing [ab|cd] ");

              computeAbcd(dLim, cLim, bLim, aLim,
                          dStrt, cStrt, bStrt, aStrt,
                          shellB, shellA,
                          aCen, bCen, p, sigmaP,
                          jMat, kMat, dMat); 
           }
        }
        return nTot;
    }

    /** MD recurrence relation steps in following two subroutines */

    private def mdRecurse(r:Vector3d, i:Int, j:Int, k:Int, m:Int) : Double {
         var res:Double = 0.0;
         val ijk = i|j|k;

         if (ijk != 0n) {
           if (ijk == 1n) {
             var m_l:Int = m;
             var r_l:Double = 1.0;
             if (i == 1n) {
                m_l++;
                r_l *= r.i;
             }
             if (j == 1n) {
                m_l++;
                r_l *= r.j;
             }
             if (k == 1n) {
               m_l++;
               r_l *= r.k;
             }
             res = r_l*zeroM(m_l);
           } else {
             if (i >= 2n) {
               res = r.i*mdRecurse(r,i-1n,j,k,m+1n)-(i-1n)*mdRecurse(r,i-2n,j,k,m+1n);
             } else if (j >= 2n) {
               res = r.j*mdRecurse(r,i,j-1n,k,m+1n)-(j-1n)*mdRecurse(r,i,j-2n,k,m+1n);
             } else if (k >= 2n) {
               res = r.k*mdRecurse(r,i,j,k-1n,m+1n)-(k-1n)*mdRecurse(r,i,j,k-2n,m+1n);
             }
           } // end if
         } else {
            res = zeroM(m);
         } // end if

         return res;
    }

    private def mdHrr(xa:Int, ya:Int, za:Int, xb:Int, yb:Int, zb:Int, xp:Int, yp:Int, zp:Int,
                      pai:Double, paj:Double, pak:Double, pbi:Double, pbj:Double, pbk:Double, sigma:Double) : Double {
        var res:Double;
        if (xp < 0n || yp < 0n || zp < 0n) {
            res = 0.0;
        } else if (xa > 0n) {
            res = mdHrr(xa-1n, ya, za, xb, yb, zb, xp-1n, yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*xp
                + mdHrr(xa-1n, ya, za, xb, yb, zb, xp   , yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*(pai)
                + mdHrr(xa-1n, ya, za, xb, yb, zb, xp+1n, yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else if (ya > 0n) {
            res = mdHrr(xa, ya-1n, za, xb, yb, zb, xp, yp-1n, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*yp
                + mdHrr(xa, ya-1n, za, xb, yb, zb, xp, yp   , zp, pai, paj, pak, pbi, pbj, pbk, sigma)*(paj)
                + mdHrr(xa, ya-1n, za, xb, yb, zb, xp, yp+1n, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else if (za > 0n) {
            res = mdHrr(xa, ya, za-1n, xb, yb, zb, xp, yp, zp-1n, pai, paj, pak, pbi, pbj, pbk, sigma)*zp
                + mdHrr(xa, ya, za-1n, xb, yb, zb, xp, yp, zp   , pai, paj, pak, pbi, pbj, pbk, sigma)*(pak)
                + mdHrr(xa, ya, za-1n, xb, yb, zb, xp, yp, zp+1n, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else if (xb > 0n) {
            res = mdHrr(xa, ya, za, xb-1n, yb, zb, xp-1n, yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*xp
                + mdHrr(xa, ya, za, xb-1n, yb, zb, xp   , yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*(pbi)
                + mdHrr(xa, ya, za, xb-1n, yb, zb, xp+1n, yp, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else if (yb > 0n) {
            res = mdHrr(xa, ya, za, xb, yb-1n, zb, xp, yp-1n, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*yp
                + mdHrr(xa, ya, za, xb, yb-1n, zb, xp, yp   , zp, pai, paj, pak, pbi, pbj, pbk, sigma)*(pbj)
                + mdHrr(xa, ya, za, xb, yb-1n, zb, xp, yp+1n, zp, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else if (zb > 0n) {
            res = mdHrr(xa, ya, za, xb, yb, zb-1n, xp, yp, zp-1n, pai, paj, pak, pbi, pbj, pbk, sigma)*zp
                + mdHrr(xa, ya, za, xb, yb, zb-1n, xp, yp, zp   , pai, paj, pak, pbi, pbj, pbk, sigma)*(pbk)
                + mdHrr(xa, ya, za, xb, yb, zb-1n, xp, yp, zp+1n, pai, paj, pak, pbi, pbj, pbk, sigma)*sigma;
        } else {
            val ptot = xp+yp+zp;
            val idx  = xp*(2n*(ptot)-xp+3n)/2n+yp;
            res = npint(ptot, idx);
        }

         return res;
    }

    /** 
     * TODO The best way to evaluate this function is to follow "GJP1991"
     * P. M. W. Gill, B. G. Johnson and J. A. Pople 
     * "Two-electron repulsion Integrals Over Gaussian s Functions"
     * International Journal of Quantum Chemistry, vol 40, 745-752 (1991)
    */    
    def computeZeroM(angMomABCD:Int, Ti:Double, Upq:Double, etai:Double) {
        // ADG1999
        var eta:Double=etai,T:Double=Ti;
        if (omega>0.0) {
            val om2=omega*omega;
            val r2=Ti/etai; // eqn 6
            eta = etai * om2 / (om2+etai); // eqn 13
            T = eta* r2; // eqn 12
        }

        // GJP1991 eqn 20-21
        if (T > TCrit) { 
            val invR2 = eta/T;
            zeroM(0) = Upq*Math.sqrt(invR2);
            for (m in 1..angMomABCD)
                zeroM(m) = zeroM(m-1) * (2*m-1) * invR2;                       
            return zeroM;
        }

        // In GJP1991, this is calculated by Modified Chebyshev Interpolation instead.
        computeGmt(angMomABCD, T);
        val twoEta = 2.0*eta;
        var etaPow : Double = Math.sqrt(twoEta);
        val UpqSQ2PI = Upq*SQ2PI;
        for(i in 0..angMomABCD) {
            zeroM(i) = UpqSQ2PI * etaPow * gmt(i);
            etaPow *= twoEta; // etaPow = twoEta^(i+0.5)
        }
        return zeroM;
    }

    private def computeGmt(angMomABCD:Int, T:Double) {
        if (T > 30.0){  // This is quivalent to GJP1991 eqn 20-21
            val invT = 1/T;
            gmt(0) = Math.sqrt(Math.PI*invT)*0.5;
           
            for (m in 1..angMomABCD) {
                gmt(m) = gmt(m-1) * (m-0.5) * invT;
            }
        } else {  
            // This is done by Modified Chebyshev Interpolation in GJP1991. 
            // However, this one is not too bad. term ~ T^n / (angMomABCD + 0.5 + n)!
            // By Sterling approximation one can see that the factorial kills term very quickly.
            // This is equivalent to eqn 29 in GJP1991. (Taylor's series)
            var denom : Double = angMomABCD + 0.5;
            var term : Double = 0.5 / denom;
            var sum : Double = term;
            while (term > THRESH)  {
                denom = (denom + 1.0);
                term  = term * T / denom;
                sum   = sum + term;
            }
            // This exp is done by Chebyshev Interpolation too in GJP1991
            val expNegT = Math.exp(-T);
            gmt(angMomABCD) = expNegT * sum;
            for (var m:Int=angMomABCD-1n; m >= 0n; m--) { // This is quivalent to GJP1991 eqn 24
                gmt(m) = (2.0 * T * gmt(m+1) + expNegT) / (2n*m + 1n);
            }
        }
    }

    private def computeRm(angMomABCD:Int, shellList:ShellList, r:Vector3d) {
        var iLim:Int = 0n;
        for(i in 0n..angMomABCD) {
            val shell = shellList.getPowers(i);
            iLim += (i+1n); // iLim = (i+1)*(i+2)/2
            for(var j:Int=0n; j<iLim; j++) {
                val powers = shell(j);
                val lp = powers.l;
                val mp = powers.m;
                val np = powers.n;
                rM(i,j) = mdRecurse(r, lp, mp, np, 0n);  // can use vrr() instead
                // Console.OUT.println(rM(i,j));
            }
        }
    }

    private def computePq(angMomAB:Int, angMomCD:Int, shellList:ShellList) {
         var pLim:Int = 0n;
         for(i in 0n..angMomAB) {
             val shellAB = shellList.getPowers(i);
             pLim += (i+1n); // pLim = (p+1)*(p+2)/2
             for (var pp:Int = 0n; pp<pLim; pp++) {
                 val powersAB = shellAB(pp);
                 val lp = powersAB.l;
                 val mp = powersAB.m;
                 val np = powersAB.n;

                 var qLim:Int = 0n;
                 for(j in 0n..angMomCD) {
                     val shellCD = shellList.getPowers(j);
                     qLim += (j+1n); // qLim = (q+1)*(q+2)/2
                     for (var qq:Int = 0n; qq<qLim; qq++) {
                         val powersCD = shellCD(qq);
                         val lq = powersCD.l;
                         val mq = powersCD.m;
                         val nq = powersCD.n;

                         val lr = lp+lq;
                         val mr = mp+mq;
                         val nr = np+nq;
                         val rtyp = lr+mr+nr;

                         val rr = lr*(2n*(lr+mr+nr)-lr+3n)/2n+mr;

                         if ((lq+mq+nq)%2n == 0n)
                            pqInts(i, pp, j, qq) =  rM(rtyp, rr);
                         else
                            pqInts(i, pp, j, qq) = -rM(rtyp, rr);
                     }
                 }
             }
         }
    }

    private def computePcd(angMomAB:Int, sigmaQ:Double, q:Point3d, dLim:Int, cLim:Int, 
                           shellD:Rail[Power], shellC:Rail[Power], dCen:Point3d, cCen:Point3d) {
         val qdi = q.i-dCen.i;
         val qdj = q.j-dCen.j;
         val qdk = q.k-dCen.k;
         val qci = q.i-cCen.i;
         val qcj = q.j-cCen.j;
         val qck = q.k-cCen.k;

         val halfSigmaQ = 0.5*sigmaQ;

         var pLim:Int = 0n;
         for(i in 0n..angMomAB) {
             pLim += (i+1n); // pLim = (p+1)*(p+2)/2
             for(pp in 0n..(pLim-1n)) {
                Array.copy(pqInts, pqInts.region.indexOf(i,pp,0,0), npint, 0L, maxam2N);

                 for(dd in 0n..(dLim-1n)) {
                     val powersD = shellD(dd);
                     val lp = powersD.l;
                     val mp = powersD.m;
                     val np = powersD.n;

                     for(cc in 0n..(cLim-1n)) {
                         val powersC = shellC(cc);
                         val lq = powersC.l;
                         val mq = powersC.m;
                         val nq = powersC.n;

                         pcdint(dd,cc,i,pp) += mdHrr(lp, mp, np, lq, mq, nq, 0n, 0n, 0n,
                                                    qdi, qdj, qdk, qci, qcj, qck, halfSigmaQ);   // can use hrr() instead

                     }
                 }
             }
         }
    }
 
    private def computeAbcd(dLim:Int, cLim:Int, bLim:Int, aLim:Int,
                            dStrt:Int, cStrt:Int, bStrt:Int, aStrt:Int,
                            shellB:Rail[Power], shellA:Rail[Power], 
                            aCen:Point3d, bCen:Point3d, p:Point3d, sigmaP:Double,
                            jMatrix:DenseMatrix, 
                            kMatrix:DenseMatrix,
                            dMatrix:DenseMatrix) {
        val pbi = p.i-bCen.i;
        val pbj = p.j-bCen.j;
        val pbk = p.k-bCen.k;

        val pai = p.i-aCen.i;
        val paj = p.j-aCen.j;
        val pak = p.k-aCen.k;

        val halfSigmaP = 0.5*sigmaP;

        for (dd in 0n..(dLim-1n)) {
            val ll = dStrt + dd;
            val normL = normFactors(ll);

            for (cc in 0n..(cLim-1n)) {
                val kk = cStrt + cc;
                val normK = normFactors(kk);

                if (kk < ll) continue;
                val kkll_st = kk*(kk+1)/2 + ll;

                Array.copy(pcdint, pcdint.region.indexOf(dd,cc,0,0), npint, 0L, maxam2N);

                for (bb in 0n..(bLim-1n)) {
                    val jj = bStrt + bb;
                    val normJ = normFactors(jj);
                    val powersB = shellB(bb);
                    val lp = powersB.l;
                    val mp = powersB.m;
                    val np = powersB.n;

                    for (aa in 0n..(aLim-1n)) {
                        val ii = aStrt + aa;

                        if (ii >= jj) {
                            val iijj_st = ii*(ii+1n)/2n + jj;

                            if (iijj_st >= kkll_st) {     
                                val normI = normFactors(ii);                              
                                val powersA = shellA(aa);
                                val lq = powersA.l;
                                val mq = powersA.m;
                                val nq = powersA.n;
                                val intVal = mdHrr(lp, mp, np, lq, mq, nq, 0n, 0n, 0n,
                                                            pbi, pbj, pbk, pai, paj, pak, halfSigmaP);  // can use hrr instead

                                val normIntVal = intVal
                                    * normL * normK
                                    * normJ * normI;
                                // printout twoeint // this is uncontacted. you have to add them algether before compare it with other programs
                                // Console.OUT.printf("2E %d %d %d %d = %e (norm = %e)\n",ii,jj,kk,ll,normIntVal,normL*normK*normJ*normI);
                                 fillJKMatrices(normIntVal, ii, jj, kk, ll, jMatrix, kMatrix, dMatrix);
                            } // end if
                        } // end if
                    } // for aa
                } // for bb
            } // for cc
        } // for dd
    }

    private @Inline def fillJKMatrices(twoEIntVal:Double,
                                ii:Int, jj:Int, kk:Int, ll:Int,
                                jMatrix:DenseMatrix, 
                                kMatrix:DenseMatrix,
                                dMatrix:DenseMatrix) {
        jMatrix(ii,jj) += dMatrix(kk,ll) * twoEIntVal;
        jMatrix(kk,ll) += dMatrix(ii,jj) * twoEIntVal;
        kMatrix(ii,kk) += dMatrix(jj,ll) * twoEIntVal;
        kMatrix(ii,ll) += dMatrix(jj,kk) * twoEIntVal;
        kMatrix(jj,kk) += dMatrix(ii,ll) * twoEIntVal;
        kMatrix(jj,ll) += dMatrix(ii,kk) * twoEIntVal;
        if (ii != jj) {
            jMatrix(jj,ii) += dMatrix(kk,ll) * twoEIntVal;
            jMatrix(kk,ll) += dMatrix(jj,ii) * twoEIntVal;
            kMatrix(jj,kk) += dMatrix(ii,ll) * twoEIntVal;
            kMatrix(jj,ll) += dMatrix(ii,kk) * twoEIntVal;
            kMatrix(ii,kk) += dMatrix(jj,ll) * twoEIntVal;
            kMatrix(ii,ll) += dMatrix(jj,kk) * twoEIntVal;
            if (kk != ll) {
                jMatrix(jj,ii) += dMatrix(ll,kk) * twoEIntVal;
                jMatrix(ll,kk) += dMatrix(jj,ii) * twoEIntVal;
                kMatrix(jj,ll) += dMatrix(ii,kk) * twoEIntVal;
                kMatrix(jj,kk) += dMatrix(ii,ll) * twoEIntVal;
                kMatrix(ii,ll) += dMatrix(jj,kk) * twoEIntVal;
                kMatrix(ii,kk) += dMatrix(jj,ll) * twoEIntVal;
            }
        }
        if (kk != ll) {
            jMatrix(ii,jj) += dMatrix(ll,kk) * twoEIntVal;
            jMatrix(ll,kk) += dMatrix(ii,jj) * twoEIntVal;
            kMatrix(ii,ll) += dMatrix(jj,kk) * twoEIntVal;
            kMatrix(ii,kk) += dMatrix(jj,ll) * twoEIntVal;
            kMatrix(jj,ll) += dMatrix(ii,kk) * twoEIntVal;
            kMatrix(jj,kk) += dMatrix(ii,ll) * twoEIntVal;
        }
        if (ii != kk || jj != ll) {
            jMatrix(kk,ll) += dMatrix(ii,jj) * twoEIntVal;
            jMatrix(ii,jj) += dMatrix(kk,ll) * twoEIntVal;
            kMatrix(kk,ii) += dMatrix(ll,jj) * twoEIntVal;
            kMatrix(kk,jj) += dMatrix(ll,ii) * twoEIntVal;
            kMatrix(ll,ii) += dMatrix(kk,jj) * twoEIntVal;
            kMatrix(ll,jj) += dMatrix(kk,ii) * twoEIntVal;
            if (ii != jj) {
                jMatrix(kk,ll) += dMatrix(jj,ii) * twoEIntVal;
                jMatrix(jj,ii) += dMatrix(kk,ll) * twoEIntVal;
                kMatrix(kk,jj) += dMatrix(ll,ii) * twoEIntVal;
                kMatrix(kk,ii) += dMatrix(ll,jj) * twoEIntVal;
                kMatrix(ll,jj) += dMatrix(kk,ii) * twoEIntVal;
                kMatrix(ll,ii) += dMatrix(kk,jj) * twoEIntVal;
                if (kk != ll) {
                    jMatrix(ll,kk) += dMatrix(jj,ii) * twoEIntVal;
                    jMatrix(jj,ii) += dMatrix(ll,kk) * twoEIntVal;
                    kMatrix(ll,jj) += dMatrix(kk,ii) * twoEIntVal;
                    kMatrix(ll,ii) += dMatrix(kk,jj) * twoEIntVal;
                    kMatrix(kk,jj) += dMatrix(ll,ii) * twoEIntVal;
                    kMatrix(kk,ii) += dMatrix(ll,jj) * twoEIntVal;
                }
            }
            if (kk != ll) {
                jMatrix(ll,kk) += dMatrix(ii,jj) * twoEIntVal;
                jMatrix(ii,jj) += dMatrix(ll,kk) * twoEIntVal;
                kMatrix(ll,ii) += dMatrix(kk,jj) * twoEIntVal;
                kMatrix(ll,jj) += dMatrix(kk,ii) * twoEIntVal;
                kMatrix(kk,ii) += dMatrix(ll,jj) * twoEIntVal;
                kMatrix(kk,jj) += dMatrix(ll,ii) * twoEIntVal;
            }
        }
    }
}

