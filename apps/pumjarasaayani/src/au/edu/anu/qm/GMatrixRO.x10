/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.util.ArrayList;
import x10.util.Team;
import x10.util.concurrent.AtomicInteger;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.vector.Point3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.util.SharedCounter;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;
import au.edu.anu.qm.ro.Integral_Pack;

/**
 * G matrix in HF calculation -- RO 
 */

public class GMatrixRO extends Matrix {

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    val roN:Int;
    val roL:Int;
    val roZ:Double;
    val nBasis:Int;
    
    val norm:Rail[Double];
    val temp:Rail[Double];
    val dk:Rail[Double];
    val munuk:Array[Double](3){rect,zeroBased};
    transient val aux:Integral_Pack;
    val jMatrix:Matrix;
    val kMatrix:Matrix;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom]) {
        super(N);
        this.bfs = bfs;
        this.mol = molecule;
        this.nBasis=N;
        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
  
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new Matrix(N);
        kMatrix = new Matrix(N);

        val roK = (roN+1)*(roL+1)*(roL+1);
        val maxam = bfs.getShellList().getMaximumAngularMomentum();
        val maxam1 = (maxam+1)*(maxam+2)/2;
        temp = new Rail[Double](maxam1*maxam1*roK);
        aux = new Integral_Pack(roN,roL);
        dk = new Rail[Double](roK); // eqn 15b in RO#7
 
        // Infinite memory code
        munuk = new Array[Double](0..(nBasis-1)*0..(nBasis-1)*0..(roK-1)); // Auxiliary integrals
    }


    public def compute(density:Density, mos:MolecularOrbitals) {
        timer.start(0);

        val N = getRowCount();

        val mosMat = mos.getMatrix();
        val denMat = density.getMatrix();
        val nOrbital = density.getNoOfOccupancies();
        val noOfAtoms = mol.getNumberOfAtoms();

        val roK = (roN+1)*(roL+1)*(roL+1);

        dk.clear();

        // Infinite memory code
        val muak = new Array[Double](0..(nBasis-1)*0..(nOrbital-1)*0..(roK-1)); // half-transformed auxiliary integrals eqn 16b in RO#7

        var mu:Int = 0; 
        var nu:Int = 0; 

        //val diagtwoe = new Matrix(N);
        //val diagtwoemat = diagtwoe.getMatrix();

        // centre a
        for(var a:Int=0; a<noOfAtoms; a++) {
            val aFunc = mol.getAtom(a).getBasisFunctions();
            val naFunc = aFunc.size();
            // basis functions on a
            for(var i:Int=0; i<naFunc; i++) {
                val iaFunc = aFunc.get(i);
                // centre b
                for(var b:Int=0; b<noOfAtoms; b++) {
                    val bFunc = mol.getAtom(b).getBasisFunctions();
                    val nbFunc = bFunc.size();
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);
                        //Console.OUT.printf("a=%d i=%d b=%d j=%d [naFunc=%d nbFunc=%d]\n", a,i,b,j,naFunc,nbFunc);
                        
                        var aaFunc:ContractedGaussian=iaFunc,bbFunc:ContractedGaussian=jbFunc;
                        val aa=iaFunc.getTotalAngularMomentum();
                        val bb=jbFunc.getTotalAngularMomentum();
                        val maxbraa = (aa+1)*(aa+2)/2; 
                        val maxbrab = (bb+1)*(bb+2)/2; 
                        // swap A and B if B has higher anglular momentum than A     
                        if (aa<bb) {
                            aaFunc=jbFunc; bbFunc=iaFunc;
                            //Console.OUT.printf("SWAP ab!\n");
                        }                  

                        // extract info from basisfunctions
                        // Note that iaFunc and jbFunc are ContractedGaussians
                        // val may not work?  must specify the same type as in cpp code?
                        val aang = aaFunc.getTotalAngularMomentum();
                        val aprimitive = aaFunc.getPrimitives();
                        val dConA = aprimitive.size;
                        val aPoint = aprimitive(0).origin;
                        val conA = new Rail[Double](dConA);
                        val zetaA = new Rail[Double](dConA);
                        for (ai in 0..(dConA-1)) {
                           conA(ai)=aprimitive(ai).coefficient;
                           zetaA(ai)=aprimitive(ai).exponent;
                           //Console.OUT.printf("a=%d i=%d ai=%d [conA(ai)=%e zetaA(ai)=%e]\n", a,i,ai,conA(ai),zetaA(ai));
                        }
                        // do the same for b

                        val bang = bbFunc.getTotalAngularMomentum();
                        val bprimitive = bbFunc.getPrimitives();
                        val dConB = bprimitive.size; 
                        val bPoint = bprimitive(0).origin; 
                        val conB = new Rail[Double](dConB); 
                        val zetaB = new Rail[Double](dConB); 
                        for (bi in 0..(dConB-1)) {
                           conB(bi)=bprimitive(bi).coefficient;
                           zetaB(bi)=bprimitive(bi).exponent;
                           //Console.OUT.printf("b=%d j=%d bi=%d [conA(bi)=%e zetaB(bi)=%e]\n", b,j,bi,conB(bi),zetaB(bi));
                        }

                        //val temp = new Rail[Double](maxbraa*maxbrab*roK); // Result for one batch
                        //Console.OUT.printf("temp size=%d ",maxbraa*maxbrab*roK);
                        //Console.OUT.printf("aang=%d bang=%d\n", aang,bang);
                        aux.genClass(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, temp);      

                        // transfer infomation from temp to munuk (Swap A and B again if necessary)
                        var ind:Int=0;
                        if (iaFunc.getTotalAngularMomentum()>=jbFunc.getTotalAngularMomentum())
                            for (var tmu:Int=mu; tmu<mu+maxbraa; tmu++) for (var tnu:Int=nu; tnu<nu+maxbrab; tnu++) for (var k:Int=0; k<roK; k++) {
                                //Console.OUT.printf("tmu=%d tnu=%d k=%d ind=%d val=%e\n",tmu,tnu,k,ind,temp(ind));
                                munuk(tmu,tnu,k)=norm(tmu)*norm(tnu)*temp(ind++);
                            }                                
                        else // Becareful... this is tricky ... maxbra are not swap 
                            for (var tnu:Int=nu; tnu<nu+maxbrab; tnu++) for (var tmu:Int=mu; tmu<mu+maxbraa; tmu++) for (var k:Int=0; k<roK; k++) {
                                //Console.OUT.printf("(Swap) tmu=%d tnu=%d k=%d ind=%d val=%e\n",tmu,tnu,k,ind,temp(ind));
                                munuk(tmu,tnu,k)=norm(tmu)*norm(tnu)*temp(ind++);
                            }

                        for (tmu in mu..(mu+maxbraa-1)) for (tnu in nu..(nu+maxbrab-1)) for (k in 0..(roK-1)) {
                           dk(k) += denMat(tmu,tnu)*munuk(tmu,tnu,k); // eqn 15b
                        }
 
                        for (tmu in mu..(mu+maxbraa-1)) for (tnu in nu..(nu+maxbrab-1)) for (aorb in 0..(nOrbital-1)) for (k in 0..(roK-1)) {
                            muak(tmu,aorb,k) += mosMat(aorb,tnu) * munuk(tmu,tnu,k); // eqn 16b the most expensive step!!!
                        }

                       // test
                       /* 
                       for (tmu in mu..(mu+maxbraa-1)) for (tnu in nu..(nu+maxbrab-1))  {
                            var intval:Double=0.;
                            for (k in 0..(roK-1)) intval+= munuk(tmu,tnu,k)*munuk(tmu,tnu,k);
                            Console.OUT.printf("mu=%d nu=%d intval=%e\n", tmu,tnu,intval);
                            diagtwoemat(tmu,tnu)=intval;
                       }*/

                        if (b!=noOfAtoms-1 || j!=nbFunc-1) nu+=maxbrab;
                        else {mu+=maxbraa; nu=0;}

                        //Console.OUT.printf("mu=%d nu=%d\n", mu,nu);

                    }
                }
            }
        }     

        val jMat = jMatrix.getMatrix();
        val kMat = kMatrix.getMatrix();
        val gMat = getMatrix();
        
        for (tmu in 0..(nBasis-1)) for (tnu in 0..(nBasis-1)) {
            var jContrib:Double = 0.0;
            var kContrib:Double = 0.0;
            for (k in 0..(roK-1))  {
                jContrib += munuk(tmu,tnu,k)*dk(k); // eqn 15a
                for (a in 0..(nOrbital-1)) {
                    kContrib += muak(tmu,a,k)*muak(tnu,a,k); // eqn16a 
                }
            }
            jMat(tmu,tnu) = jContrib;
            kMat(tmu,tnu) = kContrib;
        }
/*
        Console.OUT.println("diag2E RO");
        Console.OUT.println(diagtwoe);

        Console.OUT.println("J Mat RO");
        Console.OUT.println(jMatrix);

        Console.OUT.println("K Mat RO");
        Console.OUT.println(kMatrix);*/

        for (tmu in 0..(nBasis-1)) for (tnu in 0..(nBasis-1))
           gMat(tmu,tnu) = 2.0*jMat(tmu,tnu) - kMat(tmu,tnu); // eqn14

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }
 
}

