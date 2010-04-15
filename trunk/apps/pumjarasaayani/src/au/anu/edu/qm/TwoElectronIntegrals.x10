/**
 * TwoElectronIntegrals.x10 
 * 
 * Evaluate 2E integrals 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.ArrayList;

import x10.array.Array;

import x10x.matrix.Matrix;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;

/**
 * Old code based on original integral evaluation paper by:
 * H. Taketa, S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 *
 * New Murchie-Davidson (MD) code based on:
 * McMurchie, L. E.; Davidson, E. R. J Comput Phys, 26, 218, 1978.
 */
public class TwoElectronIntegrals { 
    global val basisFunctions:BasisFunctions!;
    global val molecule:Molecule[QMAtom]!;
    global val twoEInts:Array[Double]{self.at(this), rank==1};
    global val contractedList:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    val direct:Boolean;
    val noOfIntegrals:Int;

    global val fmt:Rail[Double]!, zeroM:Rail[Double]!;

    global val rM:Array[Double]{rank==2,self.at(this)};
    global val pqInts:Array[Double]{rank==2,self.at(this)};
    global val npint:Array[Double]{rank==2,self.at(this)};
    global val pcdint:Array[Double]{rank==3,self.at(this)};

    private val maxam:Int, maxam2:Int, maxam4:Int, maxamN:Int, maxam2M:Int, maxam2N:Int, pqdim:Int;

    public def this(basisFunctions:BasisFunctions!, molecule:Molecule[QMAtom]!, direct:Boolean) {
        this.basisFunctions = basisFunctions;
        this.molecule = molecule;
        this.direct = direct;        

        val noOfBasisFunctions = basisFunctions.getBasisFunctions().size();
        noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;

        if (!direct) {
           // allocate required memory
           twoEInts = new Array[Double]([0..noOfIntegrals]) as Array[Double]{self.at(this), rank==1};
           this.contractedList = null;
           compute2EShellPair(molecule);
        } else {
           this.contractedList = basisFunctions.getBasisFunctions();
           twoEInts = new Array[Double]([0..1]) as Array[Double]{self.at(this), rank==1};
        } // end if

    
        // scratch scpace needed for direct calculation    
        val shellList = basisFunctions.getShellList();

        // allocate scratch memory
        maxam  = shellList.getMaximumAngularMomentum();
        maxam2 = 2*maxam;
        maxam4 = 4*maxam;
        maxamN = ((maxam+1)*(maxam+2)/2);
        maxam2M  = ((maxam2+1)*(maxam2+2)/2);
        maxam2N  = ((maxam2+1)*(maxam2M+1));
        pqdim = maxam2M+1;

        // Console.OUT.println("alloc: " + maxam + " " + maxam2N);

        fmt    = Rail.make[Double](maxam4+1) as Rail[Double]!;
        zeroM  = Rail.make[Double](maxam4+1) as Rail[Double]!;

        rM     = new Array[Double]([0..maxam4+1, 0..((maxam4+1)*(maxam4+2)/2)]) as Array[Double]{rank==2,self.at(this)};
        pqInts = new Array[Double]([0..maxam2N, 0..maxam2N]) as Array[Double]{rank==2,self.at(this)};
        npint  = new Array[Double]([0..maxam2+1, 0..maxam2M+1]) as Array[Double]{rank==2,self.at(this)};
        pcdint = new Array[Double]([0..maxamN+1, 0..maxamN+1, 0..maxam2N]) as Array[Double]{rank==3,self.at(this)};

        // Console.OUT.println("alloc2: " + pcdint.region.size());
    }

    public def isDirect() = direct;
    public def getNumberOfIntegrals() : Int = noOfIntegrals;
    public def getBasisFunctions() : BasisFunctions{self.at(this)} = basisFunctions;
    public def getMolecule() = molecule;

    public def compute2E(i:Int, j:Int, k:Int, l:Int) : Double {
        return coulomb(contractedList.get(i), contractedList.get(j), 
                       contractedList.get(k), contractedList.get(l));
    }

    public def compute2E(ia:ContractedGaussian{self.at(this)}, ja:ContractedGaussian{self.at(this)}, 
                         ka:ContractedGaussian{self.at(this)}, la:ContractedGaussian{self.at(this)}) : Double {
        return coulomb(ia, ja, ka, la);
    }
 
    public interface TwoEAvailable {
        public def coulomb(i:Int, j:Int, k:Int, l:Int, twoE:Double) : void;
    }

    /** Default method: uses THO integrals. */
    public def compute2E(i:Int, j:Int, k:Int, l:Int, twoEUpdate:TwoEAvailable{self.at(this)}) : void {        
         var jij:Double = 0.0;

         val a = contractedList.get(i);
         val b = contractedList.get(j);
         val c = contractedList.get(k);
         val d = contractedList.get(l);

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

         var ii:Int, jj:Int, kk:Int, ll:Int;
         var iaExp:Double, iaCoef:Double, iaNorm:Double,
             jbExp:Double, jbCoef:Double, jbNorm:Double,
             kcExp:Double, kcCoef:Double, kcNorm:Double;

         val na = aExps.size();
         val nb = bExps.size();
         val nc = cExps.size();
         val nd = dExps.size();

         // TODO: x10 parallel
         for(ii=0; ii<na; ii++) {
             iaCoef = aCoefs.get(ii);
             iaExp  = aExps.get(ii);
             iaNorm = aNorms.get(ii);

             for(jj=0; jj<nb; jj++) {
                jbCoef = bCoefs.get(jj);
                jbExp  = bExps.get(jj);
                jbNorm = bNorms.get(jj);

                for(kk=0; kk<nc; kk++) {
                    kcCoef = cCoefs.get(kk);
                    kcExp  = cExps.get(kk);
                    kcNorm = cNorms.get(kk);

                    for(ll=0; ll<nd; ll++) {
                        jij += iaCoef * jbCoef * kcCoef
                               * dCoefs.get(ll)
                               * coulombRepulsion(
                                         aOrigin, iaNorm, aPower, iaExp,
                                         bOrigin, jbNorm, bPower, jbExp,
                                         cOrigin, kcNorm, cPower, kcExp,
                                         dOrigin,
                                         dNorms.get(ll),
                                         dPower,
                                         dExps.get(ll)
                                        ); 
                    } // end l loop
                } // end k loop
             } // end j loop
         } // end i loop

         jij = a.getNormalization() * b.getNormalization()
                 * c.getNormalization() * d.getNormalization() * jij;

         twoEUpdate.coulomb(i,j,k,l, jij);
    }

    private val sq2pi = Math.pow((2.0/Math.PI), 0.5);

    private val iidx  = Rail.make[Int](8);
    private val jjdx  = Rail.make[Int](8);
    private val kkdx  = Rail.make[Int](8);
    private val lldx  = Rail.make[Int](8);
    private val validIdx = Rail.make[Boolean](8);

    /* Note: M_D  routines mostly taken from Alistair's code, with a few changes. 
       Uses MD recurrance relations to evaluate higher angular momentum integrals.
       Direct update to GMtarix is based on the code in GMatrix.compute..() */
    public def compute2EAndRecord(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)}, 
                                  c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}, 
                                  shellList:ShellList{self.at(this)}, 
                                  jMat:Matrix{self.at(this)}, kMat:Matrix{self.at(this)},
                                  dMat:Density{self.at(this)}) : void {
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

         val jMatrix = jMat.getMatrix();
         val kMatrix = kMat.getMatrix();
         val dMatrix = dMat.getMatrix();

         val angMomAB = a.getMaximumAngularMomentum() + b.getMaximumAngularMomentum();
         val angMomCD = c.getMaximumAngularMomentum() + d.getMaximumAngularMomentum();
         val angMomABCD = angMomAB+angMomCD;

         var i:Int, j:Int, k:Int, l:Int;
         var bb:Int, aa:Int;
         var dd:Int, cc:Int;

         val radiusABSquared = a.distanceSquaredFrom(b); 
         val radiusCDSquared = c.distanceSquaredFrom(d);

         val shellA = shellList.getPowers(aAng);
         val shellB = shellList.getPowers(bAng);

         val nTot = aLim*bLim*cLim*dLim;
         val twoEInts = Rail.make[Double](nTot, (Int)=>0.0) as Rail[Double]!;

         // Console.OUT.println("New block - Allocated size: " + nInt);

         for(val aPrim in aPrims) {
           for(val bPrim in bPrims) {

             for(i=0; i<=maxamN+1; i++) 
               for(j=0; j<=maxamN+1; j++) 
                  for(k=0; k<=maxam2N; k++) 
                     pcdint(i,j,k) = 0.0;

             val aAlpha = aPrim.getExponent();
             val bAlpha = bPrim.getExponent();
             val gamma1 = aAlpha + bAlpha;
             val aCoeff = aPrim.getCoefficient();
             val bCoeff = bPrim.getCoefficient();

             // val p = gaussianProductCenter(aAlpha, aCen, bAlpha, bCen);

             val p = new Point3d(
                         (aAlpha * aCen.i + bAlpha * bCen.i) / gamma1,
                         (aAlpha * aCen.j + bAlpha * bCen.j) / gamma1,
                         (aAlpha * aCen.k + bAlpha * bCen.k) / gamma1
                       );

             val Gab = Math.exp(-aAlpha*bAlpha*radiusABSquared/gamma1);
             val Up = aCoeff*bCoeff*Gab*Math.pow((Math.PI/gamma1),1.5);
             // Console.OUT.println("Coeff: " + aCoeff + " " + bCoeff);
             // Console.OUT.println("Zeta, Gab, Up: " + gamma1 + " " + Gab + " " + Up);

             for(val cPrim in cPrims) {
               for(val dPrim in dPrims) {

                 val cAlpha = cPrim.getExponent();
                 val dAlpha = dPrim.getExponent();

                 val cCoeff = cPrim.getCoefficient();
                 val dCoeff = dPrim.getCoefficient();

                 val gamma2 = cAlpha + dAlpha;
                 val eta    = (gamma1*gamma2)/(gamma1+gamma2);

                 // val q = gaussianProductCenter(cAlpha, cCen, dAlpha, dCen);
                 val q = new Point3d(
                           (cAlpha * cCen.i + dAlpha * dCen.i) / gamma2,
                           (cAlpha * cCen.j + dAlpha * dCen.j) / gamma2,
                           (cAlpha * cCen.k + dAlpha * dCen.k) / gamma2
                         );

                 val Gcd = Math.exp(-cAlpha*dAlpha*radiusCDSquared/gamma2);
                 val Uq = cCoeff*dCoeff*Gcd*Math.pow((Math.PI/gamma2),1.5); 
                 // Console.OUT.println("Coeff: " + cCoeff + " " + dCoeff);
                 // Console.OUT.println("Zeta, Gcd, Uq: " + gamma2 + " " + Gcd + " " + Uq);

		 val xx = q.i-p.i;
		 val yy = q.j-p.j;
		 val zz = q.k-p.k;

                 val r = new Point3d(xx,yy,zz);
                 val radiusPQSquared:Double = xx*xx+yy*yy+zz*zz;  
                 val Upq = Up*Uq;
                 val T = radiusPQSquared * eta;

                 // Console.OUT.println("Computing FmT");
                 // Console.OUT.println("T value: " + T);
                 // Console.OUT.println("maxam: " + maxam);

                 // compute FmT
                 computeFmt(angMomABCD, T, fmt);
                 // computeFmtFGamma(angMomABCD, T, fmt);
        
                 // convert to GmT
                 computeGmt(angMomABCD);

                 // Console.OUT.println("Computing [0]m");

                 // Console.OUT.println("Up, Uq, Upq: " + Up + " " + Uq + " " + Upq);

                 // compute [0]m
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

        // Console.OUT.println("Filling in JK matrix");

        fillJKMatrices(dLim, cLim, bLim, aLim,
                       dStrt, cStrt, bStrt, aStrt,
                       shellList, bAng, aAng,
                       twoEInts,
                       jMatrix, kMatrix,
                       dMatrix);
    }


    // TODO: following three are duplicate methods from GMatrix, must be cleaned

    /** find unique elements and mark the onces that are not */
    private def filterUniqueElements(idx:Rail[Int]!, jdx:Rail[Int]!,
                                     kdx:Rail[Int]!, ldx:Rail[Int]!,
                                     validIdx:Rail[Boolean]!) : void {
        var i:Int, j:Int, k:Int, l:Int, m:Int, n:Int;

        for(m=0; m<8; m++) {
            i = idx(m); j = jdx(m); k = kdx(m); l = ldx(m);
            for(n=m+1; n<8; n++) {
                if (i==idx(n) && j==jdx(n) && k==kdx(n) && l==ldx(n))
                    validIdx(n) = false;
            } // end for
        } // end for
    }

    /** Set the J and K value for a given combination */
    private def setJKMatrixElements(jMatrix:Array[Double]{rank==2}, kMatrix:Array[Double]{rank==2},
                                    dMatrix:Array[Double]{rank==2},
                                    i:Int, j:Int, k:Int, l:Int, twoEIntVal:Double) : void {
        val v1 = dMatrix(k,l) * twoEIntVal;
        val v2 = dMatrix(i,j) * twoEIntVal;
        val v3 = dMatrix(j,l) * twoEIntVal;
        val v4 = dMatrix(j,k) * twoEIntVal;
        val v5 = dMatrix(i,l) * twoEIntVal;
        val v6 = dMatrix(i,k) * twoEIntVal;

        // atomic {
          jMatrix(i,j) += v1;
          jMatrix(k,l) += v2;
          kMatrix(i,k) += v3;
          kMatrix(i,l) += v4;
          kMatrix(j,k) += v5;
          kMatrix(j,l) += v6;
        // } // atomic
    }

    /** MD recurrance relation steps in following two subroutines */

    private def mdRecurse(r:Point3d, i:Int, j:Int, k:Int, m:Int) : Double {
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
         var res:Double = 0.0;

	 val al = xa|ya|za;
         val bl = xb|yb|zb;

         if (al != 0) {
           if (xa != 0 ) {
              res =   mdHrr(xa-1, ya, za, xb, yb, zb, xp-1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*xp
                        + mdHrr(xa-1, ya, za, xb, yb, zb, xp  , yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pai)
                        + mdHrr(xa-1, ya, za, xb, yb, zb, xp+1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } else if (ya != 0) {
              res =   mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp-1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*yp
                        + mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp  , zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(paj)
                        + mdHrr(xa, ya-1, za, xb, yb, zb, xp, yp+1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } else if (za != 0) {
              res =   mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp-1, pai, paj, pak, pbi, pbj, pbk, zeta2)*zp
                        + mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp  , pai, paj, pak, pbi, pbj, pbk, zeta2)*(pak)
                        + mdHrr(xa, ya, za-1, xb, yb, zb, xp, yp, zp+1, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } // end if
         } else if (bl != 0) {
           if (xb != 0 ) {
             res =   mdHrr(xa, ya, za, xb-1, yb, zb, xp-1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*xp
                     + mdHrr(xa, ya, za, xb-1, yb, zb, xp  , yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbi)
                     + mdHrr(xa, ya, za, xb-1, yb, zb, xp+1, yp, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } else if (yb != 0) {
             res =   mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp-1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*yp
                     + mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp  , zp, pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbj)
                     + mdHrr(xa, ya, za, xb, yb-1, zb, xp, yp+1, zp, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } else if (zb != 0) {
             res =   mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp-1, pai, paj, pak, pbi, pbj, pbk, zeta2)*zp
                     + mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp  , pai, paj, pak, pbi, pbj, pbk, zeta2)*(pbk)
                     + mdHrr(xa, ya, za, xb, yb, zb-1, xp, yp, zp+1, pai, paj, pak, pbi, pbj, pbk, zeta2)/(zeta2);
           } // end if
         } else if ( xp < 0 || yp < 0 || zp < 0) {
           res = 0.0; 
         } else {
           val ptot = xp+yp+zp;
           val idx  = xp*(2*(ptot)-xp+3)/2+yp;
           res  = npint(ptot, idx);
         } // end if

         return res;
    }

    private def computeGmt(angMomABCD:Int) {
         for(var i:Int=0; i<=angMomABCD; i++) {
             fmt(i) *= sq2pi;
             // Console.OUT.println(fmt(i));
         }
    }

    private def computeZeroM(angMomABCD:Int, Upq:Double, eta:Double) {
         for(var i:Int=0; i<=angMomABCD; i++) {
             zeroM(i) = Upq * Math.pow(2.0*eta, i+0.5) * fmt(i);
             // Console.OUT.println(Upq + " " + zeroM(i));
         }
    }

    private def computeRm(angMomABCD:Int, shellList:ShellList!, r:Point3d) {
         var i:Int, j:Int;

         for(i=0; i<=angMomABCD; i++) {
             val shell = shellList.getPowers(i);
             for(j=0; j<((i+1)*(i+2)/2); j++) {
                 val powers = shell(j);
                 val lp = powers.getL();
                 val mp = powers.getM();
                 val np = powers.getN();
                 rM(i,j) = mdRecurse(r, lp, mp, np, 0);  // can use vrr() instead
                 // Console.OUT.println(rM(i,j));
             }
         }
    }

    private def computePq(angMomAB:Int, angMomCD:Int, shellList:ShellList!) {
         var i:Int, j:Int, pp:Int, qq:Int;
         
         for(i=0; i<=angMomAB; i++) {
             val shellAB = shellList.getPowers(i);
             for (pp = 0; pp<((i+1)*(i+2)/2); pp++) {
                 val powersAB = shellAB(pp);
                 val lp = powersAB.getL();
                 val mp = powersAB.getM();
                 val np = powersAB.getN();

                 for(j=0; j<=angMomCD; j++) {
                     val shellCD = shellList.getPowers(j);
                     for (qq = 0; qq<((j+1)*(j+2)/2); qq++) {
                         val powersCD = shellCD(qq);
                         val lq = powersCD.getL();
                         val mq = powersCD.getM();
                         val nq = powersCD.getN();

                         val lr = lp+lq;
                         val mr = mp+mq;
                         val nr = np+nq;
                         val rtyp = lr+mr+nr;

                         val rr = lr*(2*(lr+mr+nr)-lr+3)/2+mr;

                         if ((lq+mq+nq)%2 == 0)
                            pqInts(i*pqdim+pp, j*pqdim+qq) =  rM(rtyp, rr);
                         else
                            pqInts(i*pqdim+pp, j*pqdim+qq) = -rM(rtyp, rr);

                         // Console.OUT.println(pqInts(i*pqdim+pp, j*pqdim+qq));
                     }
                 }
             }
         }
    }

    private def computePcd(angMomAB:Int, gamma2:Double, q:Point3d, dLim:Int, cLim:Int, 
                           dAng:Int, cAng:Int, dCen:Point3d, cCen:Point3d, shellList:ShellList!) {
         var i:Int, j:Int, pp:Int, k:Int, l:Int, dd:Int, cc:Int;
         val shellD = shellList.getPowers(dAng);
         val shellC = shellList.getPowers(cAng);

         val qdi = q.i-dCen.i;
         val qdj = q.j-dCen.j;
         val qdk = q.k-dCen.k;
         val qci = q.i-cCen.i;
         val qcj = q.j-cCen.j;
         val qck = q.k-cCen.k;

         for(i=0; i<=angMomAB; i++) {
             for(pp=0; pp < ((i+1)*(i+2)/2); pp++) {

                 for (k=0; k<=maxam2; k++) {
                     for (l=0; l<=maxam2M; l++) {
                         npint(k,l) = pqInts(i*pqdim+pp, k*pqdim+l);
                     }
                 }

                 for(dd = 0; dd<dLim; dd++) {
                     val powersD = shellD(dd);
                     val lp = powersD.getL();
                     val mp = powersD.getM();
                     val np = powersD.getN();

                     for(cc = 0; cc<cLim; cc++) {
                         val powersC = shellC(cc);
                         val lq = powersC.getL();
                         val mq = powersC.getM();
                         val nq = powersC.getN();

                         // Console.OUT.println("md: [" + maxam + "] " + dd + " " + cc + " " + (i*pqdim+pp));

                         pcdint(dd,cc,i*pqdim+pp) += mdHrr(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                           qdi, qdj, qdk, qci, qcj, qck, 2.0*gamma2);   // can use hrr() instead

                         // Console.OUT.println("md-done: " + dd + " " + cc + " " + i*pqdim+pp);
                     }
                 }
             }
         }
    }
 
    private def computeAbcd(dLim:Int, cLim:Int, bLim:Int, aLim:Int,
                            dStrt:Int, cStrt:Int, bStrt:Int, aStrt:Int,
                            shellList:ShellList!, bAng:Int, aAng:Int, 
                            aCen:Point3d, bCen:Point3d, p:Point3d, gamma1:Double,
                            twoEInts:Rail[Double]!) {
         var dd:Int, cc:Int, bb:Int, aa:Int, k:Int, l:Int;
         val shellB = shellList.getPowers(bAng);
         val shellA = shellList.getPowers(aAng);

         val pbi = p.i-bCen.i;
         val pbj = p.j-bCen.j;
         val pbk = p.k-bCen.k;

         val pai = p.i-aCen.i;
         val paj = p.j-aCen.j;
         val pak = p.k-aCen.k;

         var intIndx:Int = 0;

         for(dd=0; dd<dLim; dd++) {
             val ll = dStrt + dd;

             for(cc=0; cc<cLim; cc++) {
                 val kk = cStrt + cc;

                 for (k=0; k<=maxam2; k++)
                    for (l=0; l<=maxam2M; l++)
                        npint(k,l) = pcdint(dd, cc, k*pqdim+l);

                 for(bb = 0; bb<bLim; bb++) {
                     val jj = bStrt + bb;
                     val powersB = shellB(bb);
                     val lp = powersB.getL();
                     val mp = powersB.getM();
                     val np = powersB.getN();

                     for(aa = 0; aa<aLim; aa++) {
                         val ii = aStrt + aa;

                         if (ii >= jj && kk >= ll) {
                             val iijj_st = ii*(ii+1)/2 + jj;
                             val kkll_st = kk*(kk+1)/2 + ll;

                             if (iijj_st >= kkll_st) {                                   
                                 val powersA = shellA(aa);
                                 val lq = powersA.getL();
                                 val mq = powersA.getM();
                                 val nq = powersA.getN();
 
                                 twoEInts(intIndx) += mdHrr(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                            pbi, pbj, pbk, pai, paj, pak, gamma1);  // can use hrr instead
                                 intIndx++;
                             } // end if
                         } // end if
                     } // for aa
                 } // for bb
             } // for cc
         } // for d
    }

    private def fillJKMatrices(dLim:Int, cLim:Int, bLim:Int, aLim:Int,
                               dStrt:Int, cStrt:Int, bStrt:Int, aStrt:Int,
                               shellList:ShellList!, bAng:Int, aAng:Int,
                               twoEInts:Rail[Double]!,
                               jMatrix:Array[Double]{rank==2}, 
                               kMatrix:Array[Double]{rank==2},
                               dMatrix:Array[Double]{rank==2}) {
         var dd:Int, cc:Int, bb:Int, aa:Int, k:Int, l:Int;
         
         var intIndx:Int = 0;

         for(dd=0; dd<dLim; dd++) {
             val ll = dStrt + dd;
             lldx(0) = ll; lldx(1) = ll; kkdx(2) = ll; kkdx(3) = ll;
             iidx(4) = ll; iidx(5) = ll; jjdx(6) = ll; jjdx(7) = ll;

             for(cc=0; cc<cLim; cc++) {
                 val kk = cStrt + cc;
                 kkdx(0) = kk; kkdx(1) = kk; lldx(2) = kk; lldx(3) = kk;
                 jjdx(4) = kk; jjdx(5) = kk; iidx(6) = kk; iidx(7) = kk;

                 for(bb = 0; bb<bLim; bb++) {
                     val jj = bStrt + bb;

                     jjdx(0) = jj; iidx(1) = jj; iidx(2) = jj; jjdx(3) = jj;
                     lldx(4) = jj; kkdx(5) = jj; lldx(6) = jj; kkdx(7) = jj;

                     for(aa = 0; aa<aLim; aa++) {
                         val ii = aStrt + aa;
                  
                         if (ii >= jj && kk >=ll) {
                             val iijj = ii*(ii+1)/2 + jj;
                             val kkll = kk*(kk+1)/2 + ll;

                             if (iijj >= kkll) {
                                    val twoEIntVal = twoEInts(intIndx++); 
                            
                                    val v1 = dMatrix(kk,ll) * twoEIntVal;
                                    val v2 = dMatrix(ii,jj) * twoEIntVal;
                                    val v3 = dMatrix(jj,ll) * twoEIntVal;
                                    val v4 = dMatrix(jj,kk) * twoEIntVal;
                                    val v5 = dMatrix(ii,ll) * twoEIntVal;
                                    val v6 = dMatrix(ii,kk) * twoEIntVal;

                                    jMatrix(ii,jj) += v1;
                                    jMatrix(kk,ll) += v2;
                                    kMatrix(ii,kk) += v3;
                                    kMatrix(ii,ll) += v4;
                                    kMatrix(jj,kk) += v5;
                                    kMatrix(jj,ll) += v6;

                                    // special case
                                    if ((ii|jj|kk|ll) == 0) continue;
                                    iidx(0) = ii; jjdx(1) = ii; jjdx(2) = ii; iidx(3) = ii;
                                    kkdx(4) = ii; lldx(5) = ii; kkdx(6) = ii; lldx(7) = ii;

                                    // else this is symmetry unique integral, so need to
                                    // use this value for all 8 combinations
                                    // (if unique)

                                    for(var valIdx:Int=0; valIdx<8; valIdx++) validIdx(valIdx) = true;

                                    // filter unique elements
                                    filterUniqueElements(iidx, jjdx, kkdx, lldx, validIdx);

                                    // and evaluate them
                                    for(var m:Int=1; m<8; m++) {
                                        if (validIdx(m)) {
                                           // setJKMatrixElements(jMatrix, kMatrix, dMatrix, iidx(m), jjdx(m), kkdx(m), lldx(m), twoEIntVal);
                                           val ii_l = iidx(m);
                                           val jj_l = jjdx(m);
                                           val kk_l = kkdx(m);
                                           val ll_l = lldx(m);

                                           val v1_l = dMatrix(kk_l,ll_l) * twoEIntVal;
                                           val v2_l = dMatrix(ii_l,jj_l) * twoEIntVal;
                                           val v3_l = dMatrix(jj_l,ll_l) * twoEIntVal;
                                           val v4_l = dMatrix(jj_l,kk_l) * twoEIntVal;
                                           val v5_l = dMatrix(ii_l,ll_l) * twoEIntVal;
                                           val v6_l = dMatrix(ii_l,kk_l) * twoEIntVal;

                                           jMatrix(ii_l,jj_l) += v1_l;
                                           jMatrix(kk_l,ll_l) += v2_l;
                                           kMatrix(ii_l,kk_l) += v3_l;
                                           kMatrix(ii_l,ll_l) += v4_l;
                                           kMatrix(jj_l,kk_l) += v5_l;
                                           kMatrix(jj_l,ll_l) += v6_l;
                                        } // end if
                                    } // end m                    
                             } // end if
                         } // end if
                     } // for aa
                 } // for bb
             } // for cc
         } // for d
    }

    /** Method to compute all 2E integrals and store in memory for conventional method,
        not used for Direct (default) algorithm */
    protected def compute2E() : void {
        val bfs = basisFunctions.getBasisFunctions();
        val noOfBasisFunctions = bfs.size();
        

        var i:Int, j:Int, k:Int, l:Int, ij:Int, kl:Int, ijkl:Int;        
        
        var bfi:ContractedGaussian{self.at(this)}, bfj:ContractedGaussian{self.at(this)}, 
            bfk:ContractedGaussian{self.at(this)}, bfl:ContractedGaussian{self.at(this)};
        
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

    /** same as above, except loop based on atom centers */
    protected def compute2EShellPair(molecule:Molecule[QMAtom]!) : void {
        val bfs = basisFunctions.getBasisFunctions();
        val noOfBasisFunctions = bfs.size();

        val noOfAtoms = molecule.getNumberOfAtoms();
        var a:Int, b:Int, c:Int, d:Int;
        var i:Int, j:Int, k:Int, l:Int;
        var naFunc:Int, nbFunc:Int, ncFunc:Int, ndFunc:Int, twoEIndx:Int;
        var aFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}, 
            bFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}, 
            cFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}, 
            dFunc:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
        var iaFunc:ContractedGaussian{self.at(this)}, jbFunc:ContractedGaussian{self.at(this)}, 
            kcFunc:ContractedGaussian{self.at(this)}, ldFunc:ContractedGaussian{self.at(this)};

        // center a
        for(a=0; a<noOfAtoms; a++) {
            aFunc = molecule.getAtom(a).getBasisFunctions();
            naFunc = aFunc.size();
            // basis functions on a
            for(i=0; i<naFunc; i++) {
                iaFunc = aFunc.get(i);

                // center b
                for(b=0; b<=a; b++) {
                    bFunc = molecule.getAtom(b).getBasisFunctions();
                    nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(j=0; j<nbFunc; j++) {
                        jbFunc = bFunc.get(j);

                        // center c
                        for(c=0; c<noOfAtoms; c++) {
                            cFunc = molecule.getAtom(c).getBasisFunctions();
                            ncFunc = cFunc.size();
                            // basis functions on c
                            for(k=0; k<ncFunc; k++) {
                                kcFunc = cFunc.get(k);

                                // center d
                                for(d=0; d<=c; d++) {
                                    dFunc = molecule.getAtom(d).getBasisFunctions();
                                    ndFunc = (d<c) ? dFunc.size() : k+1;
                                    // basis functions on d
                                    for(l=0; l<ndFunc; l++) {
                                        ldFunc = dFunc.get(l);
                                        twoEIndx = IntegralsUtils.ijkl2intindex(
                                                           iaFunc.getIndex(),
                                                           jbFunc.getIndex(),
                                                           kcFunc.getIndex(),
                                                           ldFunc.getIndex());

                                        twoEInts(twoEIndx) = compute2E(iaFunc, jbFunc,
                                                                       kcFunc, ldFunc);
                                    } // end for (l)
                                } // end for (d)
                            } // end for (k)
                        } // end for (c)
                    }  // end for (j)
                } // end for (b)
            } // end for (i)
        } // end for (a)
    }

    public def getTwoElectronIntegrals() : Array[Double]{rank==1} = twoEInts;

    /** Return coulomb integral for a given pair of <ij|kl> contracted gaussian functions */
    private def coulomb(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)},
                        c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}) : Double {
        val la = a.getTotalAngularMomentum(), 
            lb = b.getTotalAngularMomentum(),
            lc = c.getTotalAngularMomentum(),
            ld = d.getTotalAngularMomentum();

        if (la+lb+lc+ld > 0) { return coulombFlat(a,b,c,d); }
        else                 { return coulombRec(a,b,c,d);  }
    }

    /** Return coulomb integral for a given pair of <ij|kl> contracted gaussian functions, uses non recursive routine */
    private def coulombFlat(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)},
                           c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}) : Double {
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

    private def coulombRec(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)},
                           c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}) : Double {

          val jij = contrHrr(a.getOrigin(), a.getPower(),
                             a.getCoefficients(), a.getExponents(), a.getPrimNorms(),
                             b.getOrigin(), b.getPower(),
                             b.getCoefficients(), b.getExponents(), b.getPrimNorms(),
                             c.getOrigin(), b.getPower(),
                             c.getCoefficients(), c.getExponents(), c.getPrimNorms(),
                             d.getOrigin(), d.getPower(),
                             d.getCoefficients(), d.getExponents(), d.getPrimNorms()
                            );

         return (a.getNormalization() * b.getNormalization()
                 * c.getNormalization() * d.getNormalization() * jij);
    }

    private def coulombRepulsion(
                    a:Point3d, aNorm:Double, aPower:Power, aAlpha:Double,
                    b:Point3d, bNorm:Double, bPower:Power, bAlpha:Double,
                    c:Point3d, cNorm:Double, cPower:Power, cAlpha:Double,
                    d:Point3d, dNorm:Double, dPower:Power, dAlpha:Double) : Double {

        val radiusABSquared = a.distanceSquared(b);
        val radiusCDSquared = c.distanceSquared(d);

        val p:Point3d = gaussianProductCenter(aAlpha, a, bAlpha, b);
        val q:Point3d = gaussianProductCenter(cAlpha, c, dAlpha, d);

        val radiusPQSquared = p.distanceSquared(q);

        val gamma1 = aAlpha + bAlpha;
        val gamma2 = cAlpha + dAlpha;
        val delta  = 0.25 * (1/gamma1 + 1/gamma2);

        val bx = constructBArray(
                   aPower.getL(), bPower.getL(), cPower.getL(), dPower.getL(),
                   p.i, a.i, b.i, q.i, c.i, d.i,
                   gamma1, gamma2, delta);

        val by = constructBArray(
                   aPower.getM(), bPower.getM(), cPower.getM(), dPower.getM(),
                   p.j, a.j, b.j, q.j, c.j, d.j,
                   gamma1, gamma2, delta);

        val bz = constructBArray(
                   aPower.getN(), bPower.getN(), cPower.getN(), dPower.getN(),
                   p.k, a.k, b.k, q.k, c.k, d.k,
                   gamma1, gamma2, delta);

        val nbx = bx.length;
        val nby = by.length;
        val nbz = bz.length;

        var sum:Double = 0.0;
        var i:Int, j:Int, k:Int;
        val maxam = nbx+nby+nbz;
        val fmt = Rail.make[Double](maxam+1);

        // compute FmT
        computeFmt(maxam, 0.25*radiusPQSquared/delta, fmt);

        // TODO: x10 parallel
        for(i=0; i<nbx; i++) {
            for(j=0; j<nby; j++) {
                for(k=0; k<nbz; k++) {
                    sum += bx(i) * by(j) * bz(k) * fmt(i+j+k);
                } // end for
            } // end for
        } // end for

        /**
        for(i=0; i<nbx; i++) {
            for(j=0; j<nby; j++) {
                for(k=0; k<nbz; k++) {
                    sum += bx(i) * by(j) * bz(k)
                           * IntegralsUtils.computeFGamma(
                                            i+j+k, 0.25*radiusPQSquared/delta);
                } // end for
            } // end for
        } // end for
        **/

        return (2.0 * Math.pow(Math.PI, 2.5)
                  / (gamma1 * gamma2 * Math.sqrt(gamma1+gamma2))
                  * Math.exp(-aAlpha*bAlpha*radiusABSquared/gamma1)
                  * Math.exp(-cAlpha*dAlpha*radiusCDSquared/gamma2)
                  * sum * aNorm * bNorm * cNorm * dNorm);
    }

    /** Compute the base FmT() - for evaluating integrals */
    protected def computeFmtFGamma(maxam:Int, T:Double, fmt:Rail[Double]!) {
        for(var m:Int=0; m<maxam; m++) { 
            fmt(m) = IntegralsUtils.computeFGamma(m, T);
        } // end for
    }

    /** Compute the base FmT() - for evaluating integrals */
    protected def computeFmt(maxam:Int, T:Double, fmt:Rail[Double]!) {
        var m:Int, i:Int;
        var num:Double, denom:Double, term:Double, sum:Double;
        val threshold = 1.0e-15;

        // lifted!

        if (T > 30.0){
           fmt(0) = Math.sqrt(Math.PI/T)*0.5;
           
           // TODO: x10 parallel
           for (m=1; m <= maxam; m++) {
               fmt(m) = fmt(m-1) *(m-0.5) / T;
           } // end for

        } else {
           denom = maxam + 0.5;
           term  = 0.5 / denom;
           sum   = term;
           i     = 0;
           while (term > threshold)  {
             i++;
             denom = (denom + 1.0);
             term  = term * T / denom;
             sum   = sum + term;
           } // end while

           fmt(maxam) = Math.exp(-T) * sum;
           for (m=maxam-1; m >= 0; m--) {
              fmt(m) = (2.0 * T * fmt(m+1) + Math.exp(-T)) / (2.0 * m + 1);
           } // end for
        } // end if
    }

    /**
     * recursively form the columb repulsion term using HGP, stage one: form HRR 
     * HRR (Horizontal Recurrance Relation)
     */
    protected def contrHrr(a:Point3d, aPower:Power, aCoeff:ArrayList[Double]{self.at(this)},
                           aExps:ArrayList[Double]{self.at(this)}, aNorms:ArrayList[Double]{self.at(this)},
                           b:Point3d, bPower:Power, bCoeff:ArrayList[Double]{self.at(this)},
                           bExps:ArrayList[Double]{self.at(this)}, bNorms:ArrayList[Double]{self.at(this)},
                           c:Point3d, cPower:Power, cCoeff:ArrayList[Double]{self.at(this)},
                           cExps:ArrayList[Double]{self.at(this)}, cNorms:ArrayList[Double]{self.at(this)},
                           d:Point3d, dPower:Power, dCoeff:ArrayList[Double]{self.at(this)},
                           dExps:ArrayList[Double]{self.at(this)}, dNorms:ArrayList[Double]{self.at(this)}) : Double {
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
     * VRR (Vertical Recurrance Relation) contribution
     */
    protected def contrVrr(a:Point3d, aPower:Power, aCoeff:ArrayList[Double]{self.at(this)},
                           aExps:ArrayList[Double]{self.at(this)}, aNorms:ArrayList[Double]{self.at(this)},
                           b:Point3d, bCoeff:ArrayList[Double]{self.at(this)},
                           bExps:ArrayList[Double]{self.at(this)}, bNorms:ArrayList[Double]{self.at(this)},
                           c:Point3d, cPower:Power, cCoeff:ArrayList[Double]{self.at(this)},
                           cExps:ArrayList[Double]{self.at(this)}, cNorms:ArrayList[Double]{self.at(this)},
                           d:Point3d, dCoeff:ArrayList[Double]{self.at(this)},
                           dExps:ArrayList[Double]{self.at(this)}, dNorms:ArrayList[Double]{self.at(this)}) : Double {
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

    private val sqrt2PI:Double = Math.sqrt(2.0) * Math.pow(Math.PI, 1.25); 

    /**
     * VRR (Vertical Recurrance Relation)
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
     * VRR (Vertical Recurrance Relation)
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
        val Kab  = sqrt2PI / zeta * Math.exp(-aAlpha*bAlpha / zeta*rab2);
        val rcd2 = c.distanceSquared(d);
        val Kcd  = sqrt2PI / eta * Math.exp(-cAlpha*dAlpha / eta*rcd2);
        val rpq2 = p.distanceSquared(q);
        val T    = zeta*eta / zetaPlusEta*rpq2;

        res = aNorm*bNorm*cNorm*dNorm*Kab*Kcd/Math.sqrt(zetaPlusEta)
              * IntegralsUtils.computeFGamma(m, T);
        return res;
    }

   /**
     * Construct B array.
     *
     * <i> THO eq. 2.22 </i>
     */
    private def constructBArray(l1:Int, l2:Int, l3:Int, l4:Int,
                  p:Double, a:Double, b:Double, q:Double, c:Double, d:Double,
                  g1:Double, g2:Double, delta:Double) : Rail[Double]! {
        val iMax = l1+l2+l3+l4+1;  // hold all the values (max angular momentum)
        val bArr = Rail.make[Double](iMax);

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

    /** Product of two gaussians */
    public def gaussianProductCenter(alpha1:Double, a:Point3d,  
                                    alpha2:Double, b:Point3d) : Point3d {
        val gamma:Double = alpha1 + alpha2;
        val center:Point3d =  new Point3d(
                         (alpha1 * a.i + alpha2 * b.i) / gamma,
                         (alpha1 * a.j + alpha2 * b.j) / gamma,
                         (alpha1 * a.k + alpha2 * b.k) / gamma
                       );

        return center;
    }
}

