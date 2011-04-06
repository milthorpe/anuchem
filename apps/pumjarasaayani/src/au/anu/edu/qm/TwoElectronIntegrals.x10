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
package au.anu.edu.qm;

import x10.util.ArrayList;

import x10x.matrix.Matrix;
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
 *   McMurchie, L. E.; Davidson, E. R. (1978).
 *   "One- and two-electron integrals over cartesian gaussian functions".
 *   J. Comp. Phys, 26 (2) pp. 218-231.
 */
public class TwoElectronIntegrals {
    private static val SQ2PI = Math.pow((2.0/Math.PI), 0.5); 

    private val fmt:Rail[Double], zeroM:Rail[Double];

    private val rM:Array[Double](2){rect};
    private val pqInts:Array[Double](2){rect};
    private val npint:Array[Double](2){rect};
    private val pcdint:Array[Double](3){rect};

    private val maxam:Int, maxam2:Int, maxam4:Int, maxamN:Int, maxam2M:Int, maxam2N:Int, pqdim:Int;

    /**
     * @param maxam maximum angular momentum (determines total number of integrals)
     */
    public def this(maxam : Int) {
        // allocate scratch memory
        this.maxam = maxam;
        maxam2 = 2*maxam;
        maxam4 = 4*maxam;
        maxamN = ((maxam+1)*(maxam+2)/2);
        maxam2M  = ((maxam2+1)*(maxam2+2)/2);
        maxam2N  = ((maxam2+1)*(maxam2M+1));
        pqdim = maxam2M+1;

        // Console.OUT.println("alloc: " + maxam + " " + maxam2N);

        fmt    = new Array[Double](maxam4+1);
        zeroM  = new Array[Double](maxam4+1);

        rM     = new Array[Double](0..(maxam4+1) * 0..((maxam4+1) * (maxam4+2)/2));
        pqInts = new Array[Double](0..(maxam2N)  * 0..(maxam2N));
        npint  = new Array[Double](0..(maxam2+1) * 0..(maxam2M+1));
        pcdint = new Array[Double](0..(maxamN+1) * 0..(maxamN+1) * 0..(maxam2N));

        // Console.OUT.println("alloc2: " + pcdint.region.size());
    }

    /* Note: M_D  routines mostly taken from Alistair's code, with a few changes. 
       Uses MD recurrence relations to evaluate higher angular momentum integrals.
       Direct update to GMtarix is based on the code in GMatrix.compute..() */
    public def compute2EAndRecord(a:ContractedGaussian, b:ContractedGaussian, 
                                  c:ContractedGaussian, d:ContractedGaussian, 
                                  shellList:ShellList, 
                                  jMat:Matrix, kMat:Matrix,
                                  dMat:Density) : void {
         val aPrims = a.getPrimitives();
         val bPrims = b.getPrimitives();
         val cPrims = c.getPrimitives();
         val dPrims = d.getPrimitives();

         val aCen  = a.getCenter();
         val bCen  = b.getCenter();
         val cCen  = c.getCenter();
         val dCen  = d.getCenter();

         val dAng  = d.getMaximumAngularMomentum();
         val cAng  = c.getMaximumAngularMomentum();
         val bAng  = b.getMaximumAngularMomentum();
         val aAng  = a.getMaximumAngularMomentum();

         val dLim = ((dAng+1)*(dAng+2)/2);
         val cLim = ((cAng+1)*(cAng+2)/2);
         val bLim = ((bAng+1)*(bAng+2)/2);
         val aLim = ((aAng+1)*(aAng+2)/2);

         val dStrt = d.getIntIndex();
         val cStrt = c.getIntIndex();
         val bStrt = b.getIntIndex();
         val aStrt = a.getIntIndex();

         val angMomAB = aAng + bAng;
         val angMomCD = cAng + dAng;
         val angMomABCD = angMomAB+angMomCD;

         val radiusABSquared = a.distanceSquaredFrom(b); 
         val radiusCDSquared = c.distanceSquaredFrom(d);

         val nTot = aLim*bLim*cLim*dLim;
         val twoEInts = new Array[Double](nTot);

         for([ap] in aPrims) {
            val aPrim = aPrims(ap);
           val aAlpha = aPrim.getExponent();
           val aCoeff = aPrim.getCoefficient();

           for([bp] in bPrims) {
             val bPrim = bPrims(bp);
             pcdint.fill(0.0);
             val bAlpha = bPrim.getExponent();
             val gamma1 = aAlpha + bAlpha;
             val bCoeff = bPrim.getCoefficient();

             // val p = gaussianProductCenter(aAlpha, aCen, bAlpha, bCen);

             val p = Point3d(
                         (aAlpha * aCen.i + bAlpha * bCen.i) / gamma1,
                         (aAlpha * aCen.j + bAlpha * bCen.j) / gamma1,
                         (aAlpha * aCen.k + bAlpha * bCen.k) / gamma1
                       );

             val Gab = Math.exp(-aAlpha*bAlpha*radiusABSquared/gamma1);
             val pgx = Math.PI/gamma1;
             val Up = aCoeff*bCoeff*Gab*Math.sqrt(pgx*pgx*pgx);
             // Console.OUT.println("Coeff: " + aCoeff + " " + bCoeff);
             // Console.OUT.println("Zeta, Gab, Up: " + gamma1 + " " + Gab + " " + Up);

             for([cp] in cPrims) {
                val cPrim = cPrims(cp);
               val cAlpha = cPrim.getExponent();
               val cCoeff = cPrim.getCoefficient();

               for([dp] in dPrims) {
                    val dPrim = dPrims(dp);
                 val dAlpha = dPrim.getExponent();
                 val dCoeff = dPrim.getCoefficient();

                 val gamma2 = cAlpha + dAlpha;
                 val eta    = (gamma1*gamma2)/(gamma1+gamma2);

                 // val q = gaussianProductCenter(cAlpha, cCen, dAlpha, dCen);
                 val q = Point3d(
                           (cAlpha * cCen.i + dAlpha * dCen.i) / gamma2,
                           (cAlpha * cCen.j + dAlpha * dCen.j) / gamma2,
                           (cAlpha * cCen.k + dAlpha * dCen.k) / gamma2
                         );

                 val Gcd = Math.exp(-cAlpha*dAlpha*radiusCDSquared/gamma2);
                 val qgx = Math.PI/gamma2;
                 val Uq  = cCoeff*dCoeff*Gcd*Math.sqrt(qgx*qgx*qgx); 
                 // Console.OUT.println("Coeff: " + cCoeff + " " + dCoeff);
                 // Console.OUT.println("Zeta, Gcd, Uq: " + gamma2 + " " + Gcd + " " + Uq);

                 val r = q.vector(p);
                 val radiusPQSquared = r.lengthSquared();
                 val T = radiusPQSquared * eta;

                 // Console.OUT.println("Computing FmT");
                 // Console.OUT.println("T value: " + T);
                 // Console.OUT.println("maxam: " + maxam);

                 // compute FmT
                 computeFmt(angMomABCD, T, fmt);
        
                 // convert to GmT
                 computeGmt(angMomABCD);

                 // Console.OUT.println("Computing [0]m");

                 // Console.OUT.println("Up, Uq, Upq: " + Up + " " + Uq + " " + Upq);

                 // compute [0]m
                 val Upq = Up*Uq;
                 computeZeroM(angMomABCD, Upq, eta);

                 // Console.OUT.println("Computing [r]m");
                 // Console.OUT.println("abcd_ang: " + angMomABCD);

                 // form [r]m using MD (recursion)
                 computeRm(angMomABCD, shellList, r);
 
                 // Console.OUT.println("Computing [p|q] ");

                 // form [p|q] 
                 computePq(angMomAB, angMomCD, shellList);

                 // Console.OUT.println("Computing [p|cd] ");

                 // form [p|cd]
                 computePcd(angMomAB, gamma2, q, dLim, cLim,
                            dAng, cAng, dCen, cCen, shellList);
               } // dPrim
              } // cPrim

              // form [ab|cd], 
              // Console.OUT.println("Computing [ab|cd] ");

              computeAbcd(dLim, cLim, bLim, aLim,
                          dStrt, cStrt, bStrt, aStrt,
                          shellList, bAng, aAng,
                          aCen, bCen, p, 2.0*gamma1,
                          twoEInts); 
           }
        }

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();
        val dMatrix = dMat.getMatrix();

        fillJKMatrices(dLim, cLim, bLim, aLim,
                       dStrt, cStrt, bStrt, aStrt,
                       shellList, bAng, aAng,
                       twoEInts,
                       jMatrix, kMatrix,
                       dMatrix);
    }

    /** A modification of above function */
    public def compute2EAndRecord2(a:ContractedGaussian, b:ContractedGaussian, 
                                  c:ContractedGaussian, d:ContractedGaussian, 
                                  shellList:ShellList, 
                                  jMat:Matrix, kMat:Matrix,
                                  dMat:Density,
                                  radiusABSquared:Double, 
                                  aAng:Int, bAng:Int, cAng:Int, dAng:Int, angMomAB:Int,
                                  aStrt:Int, bStrt:Int, cStrt:Int, dStrt:Int,
                                  aLim:Int, bLim:Int, abLim:Int) : void {
        val aPrims = a.getPrimitives();
        val bPrims = b.getPrimitives();
        val cPrims = c.getPrimitives();
        val dPrims = d.getPrimitives();

        val aCen  = a.getCenter();
        val bCen  = b.getCenter();
        val cCen  = c.getCenter();
        val dCen  = d.getCenter();

        val dLim = ((dAng+1)*(dAng+2)/2);
        val cLim = ((cAng+1)*(cAng+2)/2);

        val angMomCD = cAng + dAng;
        val angMomABCD = angMomAB+angMomCD;

        val radiusCDSquared = c.distanceSquaredFrom(d);

        val nTot = abLim*cLim*dLim;
        val twoEInts = new Array[Double](nTot);

        for([ap] in aPrims) {
          val aPrim = aPrims(ap);
          val aAlpha = aPrim.getExponent();
          val aCoeff = aPrim.getCoefficient();

          for([bp] in bPrims) {
             val bPrim = bPrims(bp);
             pcdint.fill(0.0);
             val bAlpha = bPrim.getExponent();
             val gamma1 = aAlpha + bAlpha;
             val bCoeff = bPrim.getCoefficient();

             // val p = gaussianProductCenter(aAlpha, aCen, bAlpha, bCen);

             val p = Point3d(
                         (aAlpha * aCen.i + bAlpha * bCen.i) / gamma1,
                         (aAlpha * aCen.j + bAlpha * bCen.j) / gamma1,
                         (aAlpha * aCen.k + bAlpha * bCen.k) / gamma1
                       );

             val Gab = Math.exp(-aAlpha*bAlpha*radiusABSquared/gamma1);
             val pgx = Math.PI/gamma1;
             val Up = aCoeff*bCoeff*Gab*Math.sqrt(pgx*pgx*pgx);
             // Console.OUT.println("Coeff: " + aCoeff + " " + bCoeff);
             // Console.OUT.println("Zeta, Gab, Up: " + gamma1 + " " + Gab + " " + Up);

             for([cp] in cPrims) {
                val cPrim = cPrims(cp);
               val cAlpha = cPrim.getExponent();
               val cCoeff = cPrim.getCoefficient();

               for([dp] in dPrims) {
                    val dPrim = dPrims(dp);
                 val dAlpha = dPrim.getExponent();
                 val dCoeff = dPrim.getCoefficient();

                 val gamma2 = cAlpha + dAlpha;
                 val eta    = (gamma1*gamma2)/(gamma1+gamma2);

                 // val q = gaussianProductCenter(cAlpha, cCen, dAlpha, dCen);
                 val q = Point3d(
                           (cAlpha * cCen.i + dAlpha * dCen.i) / gamma2,
                           (cAlpha * cCen.j + dAlpha * dCen.j) / gamma2,
                           (cAlpha * cCen.k + dAlpha * dCen.k) / gamma2
                         );

                 val Gcd = Math.exp(-cAlpha*dAlpha*radiusCDSquared/gamma2);
                 val qgx = Math.PI/gamma2;
                 val Uq  = cCoeff*dCoeff*Gcd*Math.sqrt(qgx*qgx*qgx); 
                 // Console.OUT.println("Coeff: " + cCoeff + " " + dCoeff);
                 // Console.OUT.println("Zeta, Gcd, Uq: " + gamma2 + " " + Gcd + " " + Uq);

                 val r = q.vector(p);
                 val radiusPQSquared = r.lengthSquared();  
                 val T = radiusPQSquared * eta;

                 // Console.OUT.println("Computing FmT");
                 // Console.OUT.println("T value: " + T);
                 // Console.OUT.println("maxam: " + maxam);

                 // compute FmT
                 computeFmt(angMomABCD, T, fmt);
        
                 // convert to GmT
                 computeGmt(angMomABCD);

                 // Console.OUT.println("Computing [0]m");

                 // Console.OUT.println("Up, Uq, Upq: " + Up + " " + Uq + " " + Upq);

                 // compute [0]m
                 val Upq = Up*Uq;
                 computeZeroM(angMomABCD, Upq, eta);

                 // Console.OUT.println("Computing [r]m");
                 // Console.OUT.println("abcd_ang: " + angMomABCD);

                 // form [r]m using MD (recursion)
                 computeRm(angMomABCD, shellList, r);
 
                 // Console.OUT.println("Computing [p|q] ");

                 // form [p|q] 
                 computePq(angMomAB, angMomCD, shellList);

                 // Console.OUT.println("Computing [p|cd] ");

                 // form [p|cd]
                 computePcd(angMomAB, gamma2, q, dLim, cLim,
                            dAng, cAng, dCen, cCen, shellList);
               } // dPrim
              } // cPrim

              // form [ab|cd], 
              // Console.OUT.println("Computing [ab|cd] ");

              computeAbcd(dLim, cLim, bLim, aLim,
                          dStrt, cStrt, bStrt, aStrt,
                          shellList, bAng, aAng,
                          aCen, bCen, p, 2.0*gamma1,
                          twoEInts); 
           }
        }

        val jMatrix = jMat.getMatrix();
        val kMatrix = kMat.getMatrix();
        val dMatrix = dMat.getMatrix();

        fillJKMatrices(dLim, cLim, bLim, aLim,
                       dStrt, cStrt, bStrt, aStrt,
                       shellList, bAng, aAng,
                       twoEInts,
                       jMatrix, kMatrix,
                       dMatrix);
    }

    /** MD recurrence relation steps in following two subroutines */

    private def mdRecurse(r:Vector3d, i:Int, j:Int, k:Int, m:Int) : Double {
         var res:Double = 0.0;
         val ijk = i|j|k;

         if (ijk != 0) {
           if (ijk == 1) {
             var m_l:Int = m;
             var r_l:Double = 1.0;
             if (i == 1) {
                m_l++;
                r_l *= r.i;
             }
             if (j == 1) {
                m_l++;
                r_l *= r.j;
             }
             if (k == 1) {
               m_l++;
               r_l *= r.k;
             }
             res = r_l*zeroM(m_l);
           } else {
             if (i >= 2) {
               res = r.i*mdRecurse(r,i-1,j,k,m+1)-(i-1)*mdRecurse(r,i-2,j,k,m+1);
             } else if (j >= 2 ) {
               res = r.j*mdRecurse(r,i,j-1,k,m+1)-(j-1)*mdRecurse(r,i,j-2,k,m+1);
             } else if (k >= 2 ) {
               res = r.k*mdRecurse(r,i,j,k-1,m+1)-(k-1)*mdRecurse(r,i,j,k-2,m+1);
             }
           } // end if
         } else {
            res = zeroM(m);
         } // end if

         return res;
    }

    private def mdHrr(xa:Int, ya:Int, za:Int, xb:Int, yb:Int, zb:Int, xp:Int, yp:Int, zp:Int,
                      pai:Double, paj:Double, pak:Double, pbi:Double, pbj:Double, pbk:Double, zeta2:Double) : Double {
        var res:Double;
        if (xp < 0 || yp < 0 || zp < 0) {
            res = 0.0;
        } else if (xa > 0) {
            res = mdHrr(xa-1, ya, za, xb, yb, zb, xp-1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*xp
                + mdHrr(xa-1, ya, za, xb, yb, zb, xp  , yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pai)
                + mdHrr(xa-1, ya, za, xb, yb, zb, xp+1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else if (ya > 0) {
            res = mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp-1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*yp
                + mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp  , zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(paj)
                + mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp+1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else if (za > 0) {
            res = mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp-1, pai, paj, pak, pbi, pbj, pbk, zeta2)*zp
                + mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp  , pai, paj, pak, pbi, pbj, pbk, zeta2)*(pak)
                + mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp+1, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else if (xb > 0 ) {
            res = mdHrr(xa, ya, za, xb-1, yb, zb, xp-1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*xp
                + mdHrr(xa, ya, za, xb-1, yb, zb, xp  , yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbi)
                + mdHrr(xa, ya, za, xb-1, yb, zb, xp+1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else if (yb > 0) {
            res = mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp-1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*yp
                + mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp  , zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbj)
                + mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp+1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else if (zb > 0) {
            res = mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp-1, pai, paj, pak, pbi, pbj, pbk, zeta2)*zp
                + mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp  , pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbk)
                + mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp+1, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
        } else {
            val ptot = xp+yp+zp;
            val idx  = xp*(2*(ptot)-xp+3)/2+yp;
            res = npint(ptot, idx);
        }

         return res;
    }

    private def computeGmt(angMomABCD:Int) {
         for(var i:Int=0; i<=angMomABCD; i++) {
             fmt(i) *= SQ2PI;
             // Console.OUT.println(fmt(i));
         }
    }

    private def computeZeroM(angMomABCD:Int, Upq:Double, eta:Double) {
        val twoEta = 2.0*eta;
        var etaPow : Double = Math.sqrt(twoEta);
        for(i in 0..angMomABCD) {
            zeroM(i) = Upq * etaPow * fmt(i);
            etaPow *= twoEta; // etaPow = twoEta^(i+0.5)
            // Console.OUT.println(Upq + " " + zeroM(i));
        }
    }

    private def computeRm(angMomABCD:Int, shellList:ShellList, r:Vector3d) {
         for(var i:Int=0; i<=angMomABCD; i++) {
             val shell = shellList.getPowers(i);
             val iLim  = ((i+1)*(i+2)/2);
             for(var j:Int=0; j<iLim; j++) {
                 val powers = shell(j);
                 val lp = powers.l;
                 val mp = powers.m;
                 val np = powers.n;
                 rM(i,j) = mdRecurse(r, lp, mp, np, 0);  // can use vrr() instead
                 // Console.OUT.println(rM(i,j));
             }
         }
    }

    private def computePq(angMomAB:Int, angMomCD:Int, shellList:ShellList) {
         var pLim : Int = 0;
         for(i in 0..angMomAB) {
             val shellAB = shellList.getPowers(i);
             pLim += (i+1);
             for (var pp:Int = 0; pp<pLim; pp++) {
                 val powersAB = shellAB(pp);
                 val lp = powersAB.l;
                 val mp = powersAB.m;
                 val np = powersAB.n;

                 val ipp = i*pqdim+pp;

                 var qLim : Int = 0;
                 for(j in 0..angMomCD) {
                     val shellCD = shellList.getPowers(j);
                     qLim += (j+1);
                     for (var qq:Int = 0; qq<qLim; qq++) {
                         val powersCD = shellCD(qq);
                         val lq = powersCD.l;
                         val mq = powersCD.m;
                         val nq = powersCD.n;

                         val lr = lp+lq;
                         val mr = mp+mq;
                         val nr = np+nq;
                         val rtyp = lr+mr+nr;

                         val rr = lr*(2*(lr+mr+nr)-lr+3)/2+mr;

                         val jqq = j*pqdim+qq;

                         if ((lq+mq+nq)%2 == 0)
                            pqInts(ipp, jqq) =  rM(rtyp, rr);
                         else
                            pqInts(ipp, jqq) = -rM(rtyp, rr);

                         // Console.OUT.println(pqInts(i*pqdim+pp, j*pqdim+qq));
                     }
                 }
             }
         }
    }

    private def computePcd(angMomAB:Int, gamma2:Double, q:Point3d, dLim:Int, cLim:Int, 
                           dAng:Int, cAng:Int, dCen:Point3d, cCen:Point3d, shellList:ShellList) {
         var i:Int, j:Int, pp:Int, k:Int, l:Int, dd:Int, cc:Int;
         val shellD = shellList.getPowers(dAng);
         val shellC = shellList.getPowers(cAng);

         val qdi = q.i-dCen.i;
         val qdj = q.j-dCen.j;
         val qdk = q.k-dCen.k;
         val qci = q.i-cCen.i;
         val qcj = q.j-cCen.j;
         val qck = q.k-cCen.k;

         val twoGamma = 2.0*gamma2;

         for(i=0; i<=angMomAB; i++) {
             val pLim = ((i+1)*(i+2)/2);
             for(pp=0; pp < pLim; pp++) {

                 val ipp = i*pqdim+pp;

                 for (k=0; k<=maxam2; k++) {
                     val kpq = k*pqdim;
                     for (l=0; l<=maxam2M; l++) 
                         npint(k,l) = pqInts(ipp, kpq+l);
                 }

                 for(dd = 0; dd<dLim; dd++) {
                     val powersD = shellD(dd);
                     val lp = powersD.l;
                     val mp = powersD.m;
                     val np = powersD.n;

                     for(cc = 0; cc<cLim; cc++) {
                         val powersC = shellC(cc);
                         val lq = powersC.l;
                         val mq = powersC.m;
                         val nq = powersC.n;

                         // Console.OUT.println("md: [" + maxam + "] " + dd + " " + cc + " " + (i*pqdim+pp));

                         pcdint(dd,cc,ipp) += mdHrr(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                    qdi, qdj, qdk, qci, qcj, qck, twoGamma);   // can use hrr() instead

                         // Console.OUT.println("md-done: " + dd + " " + cc + " " + i*pqdim+pp);
                     }
                 }
             }
         }
    }
 
    private def computeAbcd(dLim:Int, cLim:Int, bLim:Int, aLim:Int,
                            dStrt:Int, cStrt:Int, bStrt:Int, aStrt:Int,
                            shellList:ShellList, bAng:Int, aAng:Int, 
                            aCen:Point3d, bCen:Point3d, p:Point3d, gamma1:Double,
                            twoEInts:Array[Double](1){rect}) {
         val shellB = shellList.getPowers(bAng);
         val shellA = shellList.getPowers(aAng);

         val pbi = p.i-bCen.i;
         val pbj = p.j-bCen.j;
         val pbk = p.k-bCen.k;

         val pai = p.i-aCen.i;
         val paj = p.j-aCen.j;
         val pak = p.k-aCen.k;

         var intIndx:Int = 0;

         for (dd in 0..(dLim-1)) {
             val ll = dStrt + dd;

             for (cc in 0..(cLim-1)) {
                 val kk = cStrt + cc;

                 if (kk < ll) continue;
                 val kkll_st = kk*(kk+1)/2 + ll;

                 for (k in 0..maxam2) {
                    val kpq = k*pqdim;
                    for (l in 0..maxam2M)
                        npint(k,l) = pcdint(dd, cc, kpq+l);
                 }

                 for (bb in 0..(bLim-1)) {
                     val jj = bStrt + bb;
                     val powersB = shellB(bb);
                     val lp = powersB.l;
                     val mp = powersB.m;
                     val np = powersB.n;

                     for (aa in 0..(aLim-1)) {
                         val ii = aStrt + aa;

                         if (ii >= jj) {
                             val iijj_st = ii*(ii+1)/2 + jj;


                             if (iijj_st >= kkll_st) {                                   
                                 val powersA = shellA(aa);
                                 val lq = powersA.l;
                                 val mq = powersA.m;
                                 val nq = powersA.n;
 
                                 twoEInts(intIndx) += mdHrr(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                            pbi, pbj, pbk, pai, paj, pak, gamma1);  // can use hrr instead
                                 intIndx++;
                             } // end if
                         } // end if
                     } // for aa
                 } // for bb
             } // for cc
         } // for dd
    }

    private def fillJKMatrices(dLim:Int, cLim:Int, bLim:Int, aLim:Int,
                               dStrt:Int, cStrt:Int, bStrt:Int, aStrt:Int,
                               shellList:ShellList, bAng:Int, aAng:Int,
                               twoEInts:Array[Double](1){rect},
                               jMatrix:Array[Double](2){rect}, 
                               kMatrix:Array[Double](2){rect},
                               dMatrix:Array[Double](2){rect}) {
        // Console.OUT.println("Filling in JK matrix");
        var intIndx:Int = 0;

        for (ll in dStrt..(dStrt+dLim-1)) {
            val kStrt = Math.max(cStrt, ll);
            for (kk in kStrt..(cStrt+cLim-1)) {
                val kkll = kk*(kk+1)/2 + ll;

                for (jj in bStrt..(bStrt+bLim-1)) {
                    val iStrt = Math.max(aStrt, jj);
                    for (ii in iStrt..(aStrt+aLim-1)) {
                        val iijj = ii*(ii+1)/2 + jj;

                        if (iijj >= kkll) {
                            val twoEIntVal = twoEInts(intIndx++); 

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
                         } // end if
                     } // for aa
                 } // for bb
             } // for cc
         } // for dd

         // Console.OUT.println("\tNumber of actual integrals used in block: " + intIndx);
    }

    /** Compute the base FmT() - for evaluating integrals */
    protected def computeFmt(maxam:Int, T:Double, fmt:Rail[Double]) {
        if (T > 30.0){
            val invT = 1/T;
            fmt(0) = Math.sqrt(Math.PI*invT)*0.5;
           
            for (m in 1..maxam) {
                fmt(m) = fmt(m-1) * (m-0.5) * invT;
            }
        } else {
            var denom : Double = maxam + 0.5;
            var term : Double = 0.5 / denom;
            var sum : Double = term;
            while (term > 1.0e-15)  {
                denom = (denom + 1.0);
                term  = term * T / denom;
                sum   = sum + term;
            }

            val expNegT = Math.exp(-T);
            fmt(maxam) = expNegT * sum;
            for (var m:Int=maxam-1; m >= 0; m--) {
                fmt(m) = (2.0 * T * fmt(m+1) + expNegT) / (2*m + 1);
            }
        } // end if
    }

    /**
     * recursively form the columb repulsion term using HGP, stage one: form HRR 
     * HRR (Horizontal Recurrence Relation)
     */
    protected def contrHrr(a:Point3d, aPower:Power, aCoeff:ArrayList[Double],
                           aExps:ArrayList[Double], aNorms:ArrayList[Double],
                           b:Point3d, bPower:Power, bCoeff:ArrayList[Double],
                           bExps:ArrayList[Double], bNorms:ArrayList[Double],
                           c:Point3d, cPower:Power, cCoeff:ArrayList[Double],
                           cExps:ArrayList[Double], cNorms:ArrayList[Double],
                           d:Point3d, dPower:Power, dCoeff:ArrayList[Double],
                           dExps:ArrayList[Double], dNorms:ArrayList[Double]) : Double {
        val la = aPower.getL(), ma = aPower.getM(), na = aPower.getN();
        val lb = bPower.getL(), mb = bPower.getM(), nb = bPower.getN();
        val lc = cPower.getL(), mc = cPower.getM(), nc = cPower.getN();
        val ld = dPower.getL(), md = dPower.getM(), nd = dPower.getN();

        if (lb > 0) {
            val newBPower = Power(lb-1,mb,nb);
            return (contrHrr(a, Power(la+1,ma,na), aCoeff, aExps, aNorms, 
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.i-b.i)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (mb > 0) {
            val newBPower = Power(lb,mb-1,nb);
            return (contrHrr(a, Power(la,ma+1,na), aCoeff, aExps, aNorms,
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.j-b.j)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (nb > 0) {
            val newBPower = Power(lb,mb,nb-1);
            return (contrHrr(a, Power(la,ma,na+1), aCoeff, aExps, aNorms,
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.k-b.k)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (ld > 0) {
            val newDPower = Power(ld-1,md,nd);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms, 
                             b, bPower, bCoeff, bExps, bNorms,
                             c, Power(lc+1,mc,nc), cCoeff, cExps, cNorms,
                             d, newDPower, dCoeff, dExps, dNorms)
                   + (c.i-d.i)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms, 
                                b, bPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, newDPower, dCoeff, dExps, dNorms));
        } else if (md > 0) {
            val newDPower = Power(ld,md-1,nd);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms,
                             b, bPower, bCoeff, bExps, bNorms,
                             c, Power(lc,mc+1,nc), cCoeff, cExps, cNorms,
                             d, newDPower, dCoeff, dExps, dNorms)
                   + (c.j-d.j)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, bPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, newDPower, dCoeff, dExps, dNorms));
        } else if (nd > 0) {
            val newDPower = Power(ld,md,nd-1);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms,
                             b, bPower, bCoeff, bExps, bNorms,
                             c, Power(lc,mc,nc+1), cCoeff, cExps, cNorms,
                             d, newDPower, dCoeff, dExps, dNorms)
                + (c.k-d.k)
                    * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                               b, bPower, bCoeff, bExps, bNorms,
                               c, cPower, cCoeff, cExps, cNorms,
                               d, newDPower, dCoeff, dExps, dNorms));
        } // end if
        
        return contrVrr(a, aPower, aCoeff, aExps, aNorms,
                        b, bCoeff, bExps, bNorms,
                        c, cPower, cCoeff, cExps, cNorms,
                        d, dCoeff, dExps, dNorms);
    }

    /**
     * VRR (Vertical Recurrence Relation) contribution
     */
    protected def contrVrr(a:Point3d, aPower:Power, aCoeff:ArrayList[Double],
                           aExps:ArrayList[Double], aNorms:ArrayList[Double],
                           b:Point3d, bCoeff:ArrayList[Double],
                           bExps:ArrayList[Double], bNorms:ArrayList[Double],
                           c:Point3d, cPower:Power, cCoeff:ArrayList[Double],
                           cExps:ArrayList[Double], cNorms:ArrayList[Double],
                           d:Point3d, dCoeff:ArrayList[Double],
                           dExps:ArrayList[Double], dNorms:ArrayList[Double]) : Double {
        var res:Double = 0.0;

        var i:Int, j:Int, k:Int, l:Int;
        var iaExp:Double, iaCoef:Double, iaNorm:Double,
            jbExp:Double, jbCoef:Double, jbNorm:Double,
            kcExp:Double, kcCoef:Double, kcNorm:Double;
        
        for (i = 0; i < aExps.size(); i++) {
            iaCoef = aCoeff.get(i);
            iaExp = aExps.get(i);
            iaNorm = aNorms.get(i);

            for (j = 0; j < bExps.size(); j++) {
                jbCoef = bCoeff.get(j);
                jbExp = bExps.get(j);
                jbNorm = bNorms.get(j);

                for (k = 0; k < cExps.size(); k++) {
                    kcCoef = cCoeff.get(k);
                    kcExp = cExps.get(k);
                    kcNorm = cNorms.get(k);

                    for(l=0; l < dExps.size(); l++) {
                        res += iaCoef * jbCoef * kcCoef * dCoeff.get(l)
                                 * vrrWrapper(a, iaNorm, aPower, iaExp,
                                       b, jbNorm, jbExp,
                                       c, kcNorm, cPower, kcExp,
                                       d, dNorms.get(l), dExps.get(l), 0);
                    } // end for
                } // end for
            } // end for
        } // end for
        
        return res;
    }

    private static val SQRT2PI = Math.sqrt(2.0) * Math.pow(Math.PI, 1.25);

    /**
     * VRR (Vertical Recurrence Relation)
     */
    protected def vrrWrapper(
                         a:Point3d, aNorm:Double, aPower:Power, aAlpha:Double,
                         b:Point3d, bNorm:Double, bAlpha:Double,
                         c:Point3d, cNorm:Double, cPower:Power, cAlpha:Double,
                         d:Point3d, dNorm:Double, dAlpha:Double, m:Int) : Double {
        return vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha,
                   c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m);
    }

    /**
     * VRR (Vertical Recurrence Relation)
     */
    protected def vrr(a:Point3d, aNorm:Double, aPower:Power, aAlpha:Double,
                      b:Point3d, bNorm:Double, bAlpha:Double,
                      c:Point3d, cNorm:Double, cPower:Power, cAlpha:Double,
                      d:Point3d, dNorm:Double, dAlpha:Double, m:Int) : Double {
        var res:Double = 0.0;

        val p = gaussianProductCenter(aAlpha, a, bAlpha, b);
        val q = gaussianProductCenter(cAlpha, c, dAlpha, d);
        val zeta = aAlpha + bAlpha;
        val eta  = cAlpha + dAlpha;
        val zetaPlusEta = zeta + eta;
        val zetaByZetaPlusEta = zeta / zetaPlusEta;
        val etaByZetaPlusEta  = eta / zetaPlusEta;
        val w = gaussianProductCenter(zeta, p, eta, q);
        
        val la = aPower.getL();
        val ma = aPower.getM();
        val na = aPower.getN();
        val lc = cPower.getL();
        val mc = cPower.getM();
        val nc = cPower.getN();
        
        if (nc > 0) {
           val newCPower = Power(lc, mc, nc-1);
           res = (q.k-c.k)*vrr(a, aNorm, aPower, aAlpha,
                                         b, bNorm, bAlpha,
                                         c, cNorm, newCPower, cAlpha,
                                         d, dNorm, dAlpha, m)
               + (w.k-q.k)*vrr(a, aNorm, aPower, aAlpha,
                                         b, bNorm, bAlpha,
                                         c, cNorm, newCPower, cAlpha,
                                         d, dNorm, dAlpha, m+1);

           if (nc > 1) {
              val newCPower1 = Power(lc, mc, nc-2);
              res += 0.5*(nc-1) / eta*(vrr(a, aNorm, aPower, aAlpha,
                                           b, bNorm, bAlpha,
                                           c, cNorm, newCPower1, cAlpha,
                                           d, dNorm, dAlpha, m)
                    -zetaByZetaPlusEta*vrr(a, aNorm, aPower, aAlpha,
                                           b, bNorm, bAlpha,
                                           c, cNorm, newCPower1, cAlpha,
                                           d, dNorm, dAlpha, m+1));
           } // end if

           if (na > 0) {
              res += 0.5*na/zetaPlusEta*vrr(a, aNorm, Power(la, ma, na-1),
                                            aAlpha,
                                            b, bNorm, bAlpha,
                                            c, cNorm, newCPower,
                                            cAlpha,
                                            d, dNorm, dAlpha, m+1);
           } // end if

           return res;
        } else if (mc > 0) {
            val newCPower = Power(lc, mc-1, nc);
            res = (q.j-c.j)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m)
                + (w.j-q.j)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (mc > 1) {
               val newCPower1 = Power(lc, mc-2, nc);
               res += 0.5*(mc-1)/eta*(vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower1, cAlpha,
                                          d, dNorm, dAlpha, m)
                   -zetaByZetaPlusEta*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower1, cAlpha,
                                          d, dNorm, dAlpha, m+1));
            } // end if

            if (ma > 0) {
                res += 0.5*ma/zetaPlusEta*vrr(a, aNorm, Power(la, ma-1, na),
                                              aAlpha,
                                              b, bNorm, bAlpha,
                                              c, cNorm, newCPower,
                                              cAlpha,
                                              d, dNorm, dAlpha, m+1);
            } // end if
            
            return res;
        } else if (lc > 0) {
            val newCPower = Power(lc-1, mc, nc);
            res = (q.i-c.i)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m)
                + (w.i-q.i)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (lc > 1) {
               val newCPower1 = Power(lc-2, mc, nc);
               res += 0.5*(lc-1)/eta*(vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower1, cAlpha,
                                          d, dNorm, dAlpha, m)
                   -zetaByZetaPlusEta*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower1, cAlpha,
                                          d, dNorm, dAlpha, m+1));
            } // end if

            if (la > 0) {
                res += 0.5*la/zetaPlusEta*vrr(a, aNorm, Power(la-1, ma, na),
                                              aAlpha,
                                              b, bNorm, bAlpha,
                                              c, cNorm, newCPower,
                                              cAlpha,
                                              d, dNorm, dAlpha, m+1);
            } // end if

            return res;
        } else if (na > 0) {
            val newAPower = Power(la, ma, na-1);
            res = (p.k-a.k)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, cPower, cAlpha,
                                          d, dNorm, dAlpha, m) 
                + (w.k-p.k)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, cPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (na > 1) {
               val newAPower1 = Power(la, ma, na-2);
               res += 0.5*(na-1)/zeta*(vrr(a, aNorm, newAPower1, aAlpha,
                                           b, bNorm, bAlpha,
                                           c, cNorm, cPower, cAlpha,
                                           d, dNorm, dAlpha, m)
                     -etaByZetaPlusEta*vrr(a, aNorm, newAPower1, aAlpha,
                                           b, bNorm, bAlpha,
                                           c, cNorm, cPower, cAlpha,
                                           d, dNorm, dAlpha, m+1));
            } // end if

            return res;
        } else if (ma > 0) {
            val newAPower = Power(la, ma-1, na);
            res = (p.j-a.j)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m)
                + (w.j-p.j)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m+1);

            if (ma > 1) {
               val newAPower1 = Power(la, ma-2, na);
               res += 0.5*(ma-1)/zeta*(vrr(a, aNorm, newAPower1,
                                           aAlpha,
                                           b, aNorm, aAlpha,
                                           c, aNorm, cPower, cAlpha,
                                           d, aNorm, dAlpha, m)
                     -etaByZetaPlusEta*vrr(a, aNorm, newAPower1,
                                           aAlpha,
                                           b, aNorm, aAlpha,
                                           c, aNorm, cPower, cAlpha,
                                           d, aNorm, dAlpha, m+1));
            } // end if
            
            return res;
        } else if (la > 0) {
            val newAPower = Power(la-1, ma, na);
            res = (p.i-a.i)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m)
                + (w.i-p.i)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m+1);

            if (la > 1) {
                val newAPower1 = Power(la-2, ma, na);
                res += 0.5*(la-1)/zeta*(vrr(a, aNorm, newAPower1, aAlpha,
                                            b, aNorm, aAlpha,
                                            c, aNorm, cPower, cAlpha,
                                            d, aNorm, dAlpha, m)
                      -etaByZetaPlusEta*vrr(a, aNorm, newAPower1, aAlpha,
                                            b, aNorm, aAlpha,
                                            c, aNorm, cPower, cAlpha,
                                            d, aNorm, dAlpha, m+1));
            } // end if
            
            return res;
        } // end if

        val rab2 = a.distanceSquared(b);
        val Kab  = SQRT2PI / zeta * Math.exp(-aAlpha*bAlpha / zeta*rab2);
        val rcd2 = c.distanceSquared(d);
        val Kcd  = SQRT2PI / eta * Math.exp(-cAlpha*dAlpha / eta*rcd2);
        val rpq2 = p.distanceSquared(q);
        val T    = zeta*eta / zetaPlusEta*rpq2;

        res = aNorm*bNorm*cNorm*dNorm*Kab*Kcd/Math.sqrt(zetaPlusEta)
              * IntegralsUtils.computeFGamma(m, T);
        return res;
    }

    /** Product of two gaussians */
    private static def gaussianProductCenter(alpha1:Double, a:Point3d,  
                                    alpha2:Double, b:Point3d) : Point3d {
        val gamma:Double = alpha1 + alpha2;
        val center:Point3d =  Point3d(
                         (alpha1 * a.i + alpha2 * b.i) / gamma,
                         (alpha1 * a.j + alpha2 * b.j) / gamma,
                         (alpha1 * a.k + alpha2 * b.k) / gamma
                       );

        return center;
    }
}

