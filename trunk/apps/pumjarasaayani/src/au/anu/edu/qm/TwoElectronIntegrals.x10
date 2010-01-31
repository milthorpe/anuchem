/**
 * TwoElectronIntegrals.x10 
 * 
 * Evaluate 2E integrals 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;

public class TwoElectronIntegrals { 
    global var basisFunctions:BasisFunctions{self.at(this)};
    global var molecule:Molecule[QMAtom]{self.at(this)};
    global var twoEInts:Array[Double]{self.at(this), rank==1};
    global var contractedList:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    var direct:Boolean;
    var noOfIntegrals:Int;

    public def this() {
        // Dummy constructor
        twoEInts = Array.make[Double]([0..1]);
        contractedList = new ArrayList[ContractedGaussian{self.at(this)}]();
        basisFunctions = null;
        direct = false;
        noOfIntegrals = 0;
    }

    public def make(bfs:BasisFunctions{self.at(this)}) { 
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

    public def make(bfs:BasisFunctions{self.at(this)}, isDirect:Boolean) {
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
           compute2EShellPair();
        } else { 
           this.contractedList = basisFunctions.getBasisFunctions();
           twoEInts = Array.make[Double]([0..1]);
        } // end if
    }

    public def make(bfs:BasisFunctions{self.at(this)}, mol:Molecule[QMAtom]{self.at(this)}, isDirect:Boolean) {
        basisFunctions = bfs;
        molecule       = mol;
        direct         = isDirect;        

        val noOfBasisFunctions = basisFunctions.getBasisFunctions().size();
        noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;

        if (!direct) {
           // allocate required memory
           twoEInts = Array.make[Double]([0..noOfIntegrals]);
           this.contractedList = null;
           compute2EShellPair();
        } else {
           this.contractedList = basisFunctions.getBasisFunctions();
           twoEInts = Array.make[Double]([0..1]);
        } // end if
    }

    public def isDirect() : Boolean = direct;
    public def getNumberOfIntegrals() : Int = noOfIntegrals;
    public def getMolecule() : Molecule[QMAtom]{self.at(this)} = molecule;
    public def getBasisFunctions() : BasisFunctions{self.at(this)} = basisFunctions;

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

    /* Note: M_D  routines mostly taken from Alistair's code, with a few changes. 
       Direct update to GMtarix is based on the code in GMatrix.compute..() */
    public def compute2EAndRecord(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)}, 
                                  c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}, 
                                  idx:Int, jdx:Int, kdx:Int, ldx:Int, shellList:ShellList{self.at(this)}, 
                                  gMat:GMatrix{self.at(this)}, dMat:Density{self.at(this)}) : void {
         val aPrims = a.getPrimitives();
         val bPrims = b.getPrimitives();
         val cPrims = c.getPrimitives();
         val dPrims = d.getPrimitives();

         val gMatrix = gMat.getMatrix();
         val dMatrix = dMat.getMatrix();

         val angMomAB = a.getTotalAngularMomentum() + b.getTotalAngularMomentum();
         val angMomCD = c.getTotalAngularMomentum() + d.getTotalAngularMomentum();

         val maxam  = angMomAB+angMomCD;
         val maxam2 = 2*maxam;


         for(val aPrim in aPrims) {
           for(val bPrim in bPrims) {

             // TODO:
             val pcdint:Array[Double]{rank==3, self.at(this)} = Array.make[Double]([0..((maxam+1)*(maxam+2)/2), 0..((maxam+1)*(maxam+2)/2), 0..maxam2*((maxam2+1)*(maxam2+2)/2)]);

             for(val cPrim in cPrims) {
               for(val dPrim in dPrims) {

                 val radiusABSquared = a.distanceSquaredFrom(b);
                 val radiusCDSquared = c.distanceSquaredFrom(d);

                 val aAlpha = aPrim.getExponent();
                 val bAlpha = bPrim.getExponent();
                 val cAlpha = cPrim.getExponent();
                 val dAlpha = dPrim.getExponent();

                 val p:Point3d{self.at(this)} = gaussianProductCenter(aAlpha, a.getCenter(), bAlpha, b.getCenter());
                 val q:Point3d{self.at(this)} = gaussianProductCenter(cAlpha, c.getCenter(), dAlpha, d.getCenter());

                 val r:Point3d{self.at(this)} = new Point3d(q.sub(p));
                 val radiusPQSquared = p.distanceSquared(q);

                 val gamma1 = aAlpha + bAlpha;
                 val gamma2 = cAlpha + dAlpha;
                 val delta  = 0.25 * (1/gamma1 + 1/gamma2);   
                 val eta    = (gamma1*gamma2)/(gamma1+gamma2);  

                 val pqd    = (maxam2+1)*(maxam2+2)/2;
                 val fmt = Rail.make[Double](maxam+1);
                 val zeroM:Array[Double]{rank==1,self.at(this)} = Array.make[Double]([0..maxam+1]);
                 val rM:Array[Double]{rank==2,self.at(this)} = Array.make[Double]([0..maxam+1, 0..shellList.getNumberOfShells()]);
                 val pqInts:Array[Double]{rank==2,self.at(this)} = Array.make[Double]([0..(maxam2+1)*(pqd+1), 0..(maxam2+1)*(pqd+1)]);
                 val pqdim = maxam2+1;

                 var i:Int, j:Int, k:Int, l:Int;

                 // compute FmT
                 computeFmt(maxam, 0.25*radiusPQSquared/delta, fmt);
        
                 // convert to GmT
                 for(i=0; i<maxam; i++) fmt(i) *= sqrt2PI;

                 // compute [0]m
                 for(i=0; i<=maxam; i++)
                    zeroM(i) = radiusPQSquared * Math.pow(2.0*eta, i+0.5) * fmt(i);

                 // form [r]m using MD (recursion)
                 for(i=0; i<=maxam; i++) {
                    val shell = shellList.getShell(i).getShellPrimitives();
                    j = 0;
                    for(val shellprim in shell) {
                       val powers = (shellprim as ContractedGaussian{self.at(this)}).getPower();
                       val lp = powers.getL();
                       val mp = powers.getM();
                       val np = powers.getN();
                       rM(i,j) = mdRecurse(r, zeroM, lp, mp, np, 0);
                       j++;
                    }
                 }
 
                 // form [p|q] 
                 for(i=0; i<=angMomAB; i++) {
                    val shellAB = shellList.getShell(i).getShellPrimitives();
                    var pp:Int = 0;
                    for(val shellprimAB in shellAB) {
                      val powersAB = (shellprimAB as ContractedGaussian{self.at(this)}).getPower();
                      val lp = powersAB.getL();
                      val mp = powersAB.getM();
                      val np = powersAB.getN();
                   
                      for(j=0; j<=angMomCD; j++) {
                        val shellCD = shellList.getShell(j).getShellPrimitives();
                        var qq:Int = 0;
                        for(val shellprimCD in shellCD) {
                          val powersCD = (shellprimCD as ContractedGaussian{self.at(this)}).getPower();
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

                          qq++;
                        } 
                      } 

                      pp++;
                    }
                 }   

                 var dd:Int, cc:Int;

                 // form [p|cd]
                 for(i=0; i<=angMomAB; i++) {
                   val shellAB = shellList.getShell(i).getShellPrimitives();
                   var pp:Int = 0;
                   for(val shellprimAB in shellAB) {
                     val npint = Array.make[Double]([0..maxam2+1, 0..((maxam2+1)*(maxam2+2)/2)+1]);

                     for (k=0; k<=maxam; k++)
                       for (l=0; l<=(maxam2+1)*(maxam2+2)/2; l++)
                          npint(k,l) = pqInts(i*pqdim+pp, k*pqdim+l);

	             val dAng = dPrim.getMaximumAngularMomentum();
	             val cAng = cPrim.getMaximumAngularMomentum();

                     val shellD = shellList.getShell(dAng).getShellPrimitives();
                     dd = 0;
                     for(val shellprimD in shellD) {
                        val powersD = (shellprimD as ContractedGaussian{self.at(this)}).getPower();
                        val lp = powersD.getL();
                        val mp = powersD.getM();
                        val np = powersD.getN();

                        val shellC = shellList.getShell(cAng).getShellPrimitives();
                        cc = 0;
                        for(val shellprimC in shellC) {
                           val powersC = (shellprimC as ContractedGaussian{self.at(this)}).getPower();
                           val lq = powersC.getL();
                           val mq = powersC.getM();
                           val nq = powersC.getN();

                           pcdint(dd,cc,i*pqdim+pp) += mdr1(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                          d.getCenter(), c.getCenter(),
                                                          q, eta, npint);
                           cc++;
                        }

                        dd++;
                     }

                     pp++;
                   }
                 }

                 // TODO:
                 // form [ab|cd], and update the GMatrix elements accordingly
                 val dAng = dPrim.getMaximumAngularMomentum();
                 val cAng = cPrim.getMaximumAngularMomentum();
                 val bAng = bPrim.getMaximumAngularMomentum();
                 val aAng = aPrim.getMaximumAngularMomentum();

                 var bb:Int, aa:Int;

                 for(dd=0; dd<((dAng+1)*(dAng+2)/2); dd++) {
                     for(cc=0; cc<((cAng+1)*(cAng+2)/2); cc++) {
                        val npint = Array.make[Double]([0..maxam2+1, 0..((maxam2+1)*(maxam2+2)/2)+1]);

                        for (k=0; k<=maxam; k++)
                          for (l=0; l<=(maxam2+1)*(maxam2+2)/2; l++)
                            npint(k,l) = pcdint(dd, cc, k*pqdim+l);

                        val shellB = shellList.getShell(bAng).getShellPrimitives();
                        bb = 0;
                        for(val shellprimB in shellB) {
                          val powersB = (shellprimB as ContractedGaussian{self.at(this)}).getPower();
                          val lp = powersB.getL();
                          val mp = powersB.getM();
                          val np = powersB.getN();

                          val shellA = shellList.getShell(aAng).getShellPrimitives();
                          aa = 0;
                          for(val shellprimA in shellA) {
                             val powersA = (shellprimA as ContractedGaussian{self.at(this)}).getPower();
                             val lq = powersA.getL();
                             val mq = powersA.getM();
                             val nq = powersA.getN();
 
                             val ii = a.getIntIndex() + aa;
                             val jj = b.getIntIndex() + bb;
                             val kk = c.getIntIndex() + cc;
                             val ll = d.getIntIndex() + dd;

                             if (ii >= jj && kk >=ll) {
                                val iijj = ii*(ii+1)/2 + jj;
                                val kkll = kk*(kk+1)/2 + ll;

                                if (iijj >= kkll) {
                                    val twoEIntVal = mdr1(lp, mp, np, lq, mq, nq, 0, 0, 0,
                                                       b.getCenter(), a.getCenter(),
                                                       p, eta, npint);
                                    val twoEIntVal2    = twoEIntVal + twoEIntVal;
                                    val twoEIntValHalf = 0.5 * twoEIntVal;

                                    val v1 = dMatrix(kk,ll) * twoEIntVal2;
                                    val v2 = dMatrix(ii,jj) * twoEIntVal2;
                                    val v3 = dMatrix(jj,ll) * twoEIntValHalf;
                                    val v4 = dMatrix(jj,kk) * twoEIntValHalf;
                                    val v5 = dMatrix(ii,ll) * twoEIntValHalf;
                                    val v6 = dMatrix(ii,kk) * twoEIntValHalf;

                                    // atomic {
                                      gMatrix(ii,jj) += v1;
                                      gMatrix(kk,ll) += v2;
                                      gMatrix(ii,kk) -= v3;
                                      gMatrix(ii,ll) -= v4;
                                      gMatrix(jj,kk) -= v5;
                                      gMatrix(jj,ll) -= v6;
                                    // } // atomic
                                }
                             }

                             aa++;
                          }

                          bb++;
                        }

                     }
                 }
              }
            }
          }
        }
    }

    protected def mdRecurse(r:Point3d{self.at(this)}, 
                            zeroM:Array[Double]{rank==1,self.at(this)},
                            i:Int, j:Int, k:Int, m:Int) : Double {
         var res:Double;

         if (i >= 2) {
           res = r.i*mdRecurse(r,zeroM,i-1,j,k,m+1)-(i-1)*mdRecurse(r,zeroM,i-2,j,k,m+1);
         } else if (j >= 2 ) {
           res = r.j*mdRecurse(r,zeroM,i,j-1,k,m+1)-(j-1)*mdRecurse(r,zeroM,i,j-2,k,m+1);
         } else if (k >= 2 ) {
           res = r.k*mdRecurse(r,zeroM,i,j,k-1,m+1)-(k-1)*mdRecurse(r,zeroM,i,j,k-2,m+1);
         } else if (i == 1 ) {
           res = r.i*mdRecurse(r,zeroM,i-1,j,k,m+1);
         } else if (j == 1 ) {
           res = r.j*mdRecurse(r,zeroM,i,j-1,k,m+1);
         } else if (k == 1 ) {
           res = r.k*mdRecurse(r,zeroM,i,j,k-1,m+1);
         } else {
           res = zeroM(m);
         } // end if

         return res;
    }

    protected def mdr1(xa:Int, ya:Int, za:Int, xb:Int, yb:Int, zb:Int, xp:Int, yp:Int, zp:Int,
                          coorda:Point3d{self.at(this)}, coordb:Point3d{self.at(this)}, coordp:Point3d{self.at(this)}, zeta:Double,
                          pint:Array[Double]{rank==2, self.at(this)}) : Double {
         var res:Double;
         var ptot:Int, idx:Int;

         if (xa != 0 ){
           res =   mdr1(xa-1, ya, za, xb, yb, zb, xp-1, yp, zp, coorda, coordb, coordp, zeta, pint)*xp
                    + mdr1(xa-1, ya, za, xb, yb, zb, xp  , yp, zp, coorda, coordb, coordp, zeta, pint)*(coordp.i-coorda.i)
                    + mdr1(xa-1, ya, za, xb, yb, zb, xp+1, yp, zp, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if (ya != 0) {
           res =   mdr1(xa, ya-1, za, xb, yb, zb, xp, yp-1, zp, coorda, coordb, coordp, zeta, pint)*yp
                    + mdr1(xa, ya-1, za, xb, yb, zb, xp, yp  , zp, coorda, coordb, coordp, zeta, pint)*(coordp.j-coorda.j)
                    + mdr1(xa, ya-1, za, xb, yb, zb, xp, yp+1, zp, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if (za != 0) {
           res =   mdr1(xa, ya, za-1, xb, yb, zb, xp, yp, zp-1, coorda, coordb, coordp, zeta, pint)*zp
                    + mdr1(xa, ya, za-1, xb, yb, zb, xp, yp, zp  , coorda, coordb, coordp, zeta, pint)*(coordp.k-coorda.k)
                    + mdr1(xa, ya, za-1, xb, yb, zb, xp, yp, zp+1, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if (xb != 0 ) {
           res =   mdr1(xa, ya, za, xb-1, yb, zb, xp-1, yp, zp, coorda, coordb, coordp, zeta, pint)*xp
                    + mdr1(xa, ya, za, xb-1, yb, zb, xp  , yp, zp, coorda, coordb, coordp, zeta, pint)*(coordp.i-coordb.i)
                    + mdr1(xa, ya, za, xb-1, yb, zb, xp+1, yp, zp, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if (yb != 0) {
           res =   mdr1(xa, ya, za, xb, yb-1, zb, xp, yp-1, zp, coorda, coordb, coordp, zeta, pint)*yp
                    + mdr1(xa, ya, za, xb, yb-1, zb, xp, yp  , zp, coorda, coordb, coordp, zeta, pint)*(coordp.j-coordb.j)
                    + mdr1(xa, ya, za, xb, yb-1, zb, xp, yp+1, zp, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if (zb != 0) {
           res =   mdr1(xa, ya, za, xb, yb, zb-1, xp, yp, zp-1, coorda, coordb, coordp, zeta, pint)*zp
                    + mdr1(xa, ya, za, xb, yb, zb-1, xp, yp, zp  , coorda, coordb, coordp, zeta, pint)*(coordp.k-coordb.k)
                    + mdr1(xa, ya, za, xb, yb, zb-1, xp, yp, zp+1, coorda, coordb, coordp, zeta, pint)/(2.0*zeta);
         } else if ( xp < 0 || yp < 0 || zp < 0) {
           res = 1.0;
         } else {
           ptot=xp+yp+zp;
           idx = xp*(2*(xp+yp+zp)-xp+3)/2+yp;
           res = pint(ptot, idx);
         } // end if

         return res;
    }


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

    protected def compute2EShellPair() : void {
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


    private def coulomb(a:ContractedGaussian{self.at(this)}, b:ContractedGaussian{self.at(this)},
                        c:ContractedGaussian{self.at(this)}, d:ContractedGaussian{self.at(this)}) : Double {
        val la = a.getTotalAngularMomentum(), 
            lb = b.getTotalAngularMomentum(),
            lc = c.getTotalAngularMomentum(),
            ld = d.getTotalAngularMomentum();

        if (la+lb+lc+ld > 0) { return coulombFlat(a,b,c,d); }
        else                 { return coulombRec(a,b,c,d);  }
    }

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

          val jij = contrHrr(a.getOrigin() as Point3d{self.at(this)}, a.getPower() as Power{self.at(this)},
                             a.getCoefficients(), a.getExponents(), a.getPrimNorms(),
                             b.getOrigin() as Point3d{self.at(this)}, b.getPower() as Power{self.at(this)},
                             b.getCoefficients(), b.getExponents(), b.getPrimNorms(),
                             c.getOrigin() as Point3d{self.at(this)}, b.getPower() as Power{self.at(this)},
                             c.getCoefficients(), c.getExponents(), c.getPrimNorms(),
                             d.getOrigin() as Point3d{self.at(this)}, d.getPower() as Power{self.at(this)},
                             d.getCoefficients(), d.getExponents(), d.getPrimNorms()
                            );

         return (a.getNormalization() * b.getNormalization()
                 * c.getNormalization() * d.getNormalization() * jij);
    }

    private def coulombRepulsion(
                    a:Point3d!, aNorm:Double, aPower:Power!, aAlpha:Double,
                    b:Point3d!, bNorm:Double, bPower:Power!, bAlpha:Double,
                    c:Point3d!, cNorm:Double, cPower:Power!, cAlpha:Double,
                    d:Point3d!, dNorm:Double, dPower:Power!, dAlpha:Double) : Double {

        val radiusABSquared = a.distanceSquared(b);
        val radiusCDSquared = c.distanceSquared(d);

        val p:Point3d{self.at(this)} = gaussianProductCenter(aAlpha, a, bAlpha, b);
        val q:Point3d{self.at(this)} = gaussianProductCenter(cAlpha, c, dAlpha, d);

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

        var sum:Double = 0.0;
        var i:Int, j:Int, k:Int;
        val nbx = bx.length;
        val nby = by.length;
        val nbz = bz.length;
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
    protected def contrHrr(a:Point3d{self.at(this)}, aPower:Power{self.at(this)}, aCoeff:ArrayList[Double]{self.at(this)},
                           aExps:ArrayList[Double]{self.at(this)}, aNorms:ArrayList[Double]{self.at(this)},
                           b:Point3d{self.at(this)}, bPower:Power{self.at(this)}, bCoeff:ArrayList[Double]{self.at(this)},
                           bExps:ArrayList[Double]{self.at(this)}, bNorms:ArrayList[Double]{self.at(this)},
                           c:Point3d{self.at(this)}, cPower:Power{self.at(this)}, cCoeff:ArrayList[Double]{self.at(this)},
                           cExps:ArrayList[Double]{self.at(this)}, cNorms:ArrayList[Double]{self.at(this)},
                           d:Point3d{self.at(this)}, dPower:Power{self.at(this)}, dCoeff:ArrayList[Double]{self.at(this)},
                           dExps:ArrayList[Double]{self.at(this)}, dNorms:ArrayList[Double]{self.at(this)}) : Double {
        val la = aPower.getL(), ma = aPower.getM(), na = aPower.getN();
        val lb = bPower.getL(), mb = bPower.getM(), nb = bPower.getN();
        val lc = cPower.getL(), mc = cPower.getM(), nc = cPower.getN();
        val ld = dPower.getL(), md = dPower.getM(), nd = dPower.getN();

        if (lb > 0) {
            val newBPower = new Power(lb-1,mb,nb);
            return (contrHrr(a, new Power(la+1,ma,na), aCoeff, aExps, aNorms, 
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.i-b.i)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (mb > 0) {
            val newBPower = new Power(lb,mb-1,nb);
            return (contrHrr(a, new Power(la,ma+1,na), aCoeff, aExps, aNorms,
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.j-b.j)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (nb > 0) {
            val newBPower = new Power(lb,mb,nb-1);
            return (contrHrr(a, new Power(la,ma,na+1), aCoeff, aExps, aNorms,
                             b, newBPower, bCoeff, bExps, bNorms,
                             c, cPower, cCoeff, cExps, cNorms,
                             d, dPower, dCoeff, dExps, dNorms)
                   + (a.k-b.k)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, newBPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, dPower, dCoeff, dExps, dNorms));
        } else if (ld > 0) {
            val newDPower = new Power(ld-1,md,nd);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms, 
                             b, bPower, bCoeff, bExps, bNorms,
                             c, new Power(lc+1,mc,nc), cCoeff, cExps, cNorms,
                             d, newDPower, dCoeff, dExps, dNorms)
                   + (c.i-d.i)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms, 
                                b, bPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, newDPower, dCoeff, dExps, dNorms));
        } else if (md > 0) {
            val newDPower = new Power(ld,md-1,nd);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms,
                             b, bPower, bCoeff, bExps, bNorms,
                             c, new Power(lc,mc+1,nc), cCoeff, cExps, cNorms,
                             d, newDPower, dCoeff, dExps, dNorms)
                   + (c.j-d.j)
                     * contrHrr(a, aPower, aCoeff, aExps, aNorms,
                                b, bPower, bCoeff, bExps, bNorms,
                                c, cPower, cCoeff, cExps, cNorms,
                                d, newDPower, dCoeff, dExps, dNorms));
        } else if (nd > 0) {
            val newDPower = new Power(ld,md,nd-1);
            return (contrHrr(a, aPower, aCoeff, aExps, aNorms,
                             b, bPower, bCoeff, bExps, bNorms,
                             c, new Power(lc,mc,nc+1), cCoeff, cExps, cNorms,
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
    protected def contrVrr(a:Point3d{self.at(this)}, aPower:Power{self.at(this)}, aCoeff:ArrayList[Double]{self.at(this)},
                           aExps:ArrayList[Double]{self.at(this)}, aNorms:ArrayList[Double]{self.at(this)},
                           b:Point3d{self.at(this)}, bCoeff:ArrayList[Double]{self.at(this)},
                           bExps:ArrayList[Double]{self.at(this)}, bNorms:ArrayList[Double]{self.at(this)},
                           c:Point3d{self.at(this)}, cPower:Power{self.at(this)}, cCoeff:ArrayList[Double]{self.at(this)},
                           cExps:ArrayList[Double]{self.at(this)}, cNorms:ArrayList[Double]{self.at(this)},
                           d:Point3d{self.at(this)}, dCoeff:ArrayList[Double]{self.at(this)},
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
                         a:Point3d{self.at(this)}, aNorm:Double, aPower:Power{self.at(this)}, aAlpha:Double,
                         b:Point3d{self.at(this)}, bNorm:Double, bAlpha:Double,
                         c:Point3d{self.at(this)}, cNorm:Double, cPower:Power{self.at(this)}, cAlpha:Double,
                         d:Point3d{self.at(this)}, dNorm:Double, dAlpha:Double, m:Int) : Double {
        return vrr(a, aNorm, aPower, aAlpha, b, bNorm, bAlpha,
                   c, cNorm, cPower, cAlpha, d, dNorm, dAlpha, m);
    }

    /**
     * VRR (Vertical Recurrance Relation)
     */
    protected def vrr(a:Point3d{self.at(this)}, aNorm:Double, aPower:Power{self.at(this)}, aAlpha:Double,
                      b:Point3d{self.at(this)}, bNorm:Double, bAlpha:Double,
                      c:Point3d{self.at(this)}, cNorm:Double, cPower:Power{self.at(this)}, cAlpha:Double,
                      d:Point3d{self.at(this)}, dNorm:Double, dAlpha:Double, m:Int) : Double {
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
           val newCPower = new Power(lc, mc, nc-1);
           res = (q.k-c.k)*vrr(a, aNorm, aPower, aAlpha,
                                         b, bNorm, bAlpha,
                                         c, cNorm, newCPower, cAlpha,
                                         d, dNorm, dAlpha, m)
               + (w.k-q.k)*vrr(a, aNorm, aPower, aAlpha,
                                         b, bNorm, bAlpha,
                                         c, cNorm, newCPower, cAlpha,
                                         d, dNorm, dAlpha, m+1);

           if (nc > 1) {
              val newCPower1 = new Power(lc, mc, nc-2);
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
              res += 0.5*na/zetaPlusEta*vrr(a, aNorm, new Power(la, ma, na-1),
                                            aAlpha,
                                            b, bNorm, bAlpha,
                                            c, cNorm, newCPower,
                                            cAlpha,
                                            d, dNorm, dAlpha, m+1);
           } // end if

           return res;
        } else if (mc > 0) {
            val newCPower = new Power(lc, mc-1, nc);
            res = (q.j-c.j)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m)
                + (w.j-q.j)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (mc > 1) {
               val newCPower1 = new Power(lc, mc-2, nc);
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
                res += 0.5*ma/zetaPlusEta*vrr(a, aNorm, new Power(la, ma-1, na),
                                              aAlpha,
                                              b, bNorm, bAlpha,
                                              c, cNorm, newCPower,
                                              cAlpha,
                                              d, dNorm, dAlpha, m+1);
            } // end if
            
            return res;
        } else if (lc > 0) {
            val newCPower = new Power(lc-1, mc, nc);
            res = (q.i-c.i)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m)
                + (w.i-q.i)*vrr(a, aNorm, aPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, newCPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (lc > 1) {
               val newCPower1 = new Power(lc-2, mc, nc);
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
                res += 0.5*la/zetaPlusEta*vrr(a, aNorm, new Power(la-1, ma, na),
                                              aAlpha,
                                              b, bNorm, bAlpha,
                                              c, cNorm, newCPower,
                                              cAlpha,
                                              d, dNorm, dAlpha, m+1);
            } // end if

            return res;
        } else if (na > 0) {
            val newAPower = new Power(la, ma, na-1);
            res = (p.k-a.k)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, cPower, cAlpha,
                                          d, dNorm, dAlpha, m) 
                + (w.k-p.k)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, bNorm, bAlpha,
                                          c, cNorm, cPower, cAlpha,
                                          d, dNorm, dAlpha, m+1);

            if (na > 1) {
               val newAPower1 = new Power(la, ma, na-2);
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
            val newAPower = new Power(la, ma-1, na);
            res = (p.j-a.j)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m)
                + (w.j-p.j)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m+1);

            if (ma > 1) {
               val newAPower1 = new Power(la, ma-2, na);
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
            val newAPower = new Power(la-1, ma, na);
            res = (p.i-a.i)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m)
                + (w.i-p.i)*vrr(a, aNorm, newAPower, aAlpha,
                                          b, aNorm, aAlpha,
                                          c, aNorm, cPower, cAlpha,
                                          d, aNorm, dAlpha, m+1);

            if (la > 1) {
                val newAPower1 = new Power(la-2, ma, na);
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
    public def gaussianProductCenter(alpha1:Double, a:Point3d{self.at(this)},  
                                    alpha2:Double, b:Point3d{self.at(this)}) : Point3d{self.at(this)} {
        val gamma:Double = alpha1 + alpha2;
        val center:Point3d{self.at(this)} =  new Point3d(
                         (alpha1 * a.i + alpha2 * b.i) / gamma,
                         (alpha1 * a.j + alpha2 * b.j) / gamma,
                         (alpha1 * a.k + alpha2 * b.k) / gamma
                       );

        return center;
    }
}

