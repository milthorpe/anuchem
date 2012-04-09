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

import x10.matrix.DenseMatrix;

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

public class GMatrixRO extends DenseMatrix {

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
    val jMatrix:DenseMatrix;
    val kMatrix:DenseMatrix;

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom]) {
        super(N,N);
        this.bfs = bfs;
        this.mol = molecule;
        this.nBasis=N;
        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        this.roL=jd.roL;
        this.roZ=jd.roZ;
  
        this.norm = bfs.getNormalizationFactors();
        jMatrix = new DenseMatrix(N,N);
        kMatrix = new DenseMatrix(N,N);

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

        val nOrbital = density.getNoOfOccupancies();
        val noOfAtoms = mol.getNumberOfAtoms();

        val roK = (roN+1)*(roL+1)*(roL+1);

        dk.clear();

        // Infinite memory code
        val muak = new Array[Double](0..(nBasis-1)*0..(nOrbital-1)*0..(roK-1)); // half-transformed auxiliary integrals eqn 16b in RO#7

        var mu:Int = 0; 
        var nu:Int = 0; 

        //val diagtwoe = new DenseMatrix(N);
        //val diagtwoemat = diagtwoe.DenseMatrix();

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
                        val aPoint = aaFunc.origin;
                        val zetaA = aaFunc.exponents;
                        val conA = aaFunc.coefficients;
                        val dConA = conA.size;

                        val bang = bbFunc.getTotalAngularMomentum();
                        val bPoint = bbFunc.origin; 
                        val zetaB = bbFunc.exponents;                        
                        val conB = bbFunc.coefficients; 
                        val dConB = conB.size;

                        //Console.OUT.printf("aang=%d bang=%d\n", aang,bang);
                        aux.genClass(aang, bang, aPoint, bPoint, zetaA, zetaB, conA, conB, dConA, dConB, temp);      

                        // transfer infomation from temp to munuk (Swap A and B again if necessary) and normalise
                        var ind:Int=0;
                        if (iaFunc.getTotalAngularMomentum()>=jbFunc.getTotalAngularMomentum())
                            for (var tmu:Int=mu; tmu<mu+maxbraa; tmu++) for (var tnu:Int=nu; tnu<nu+maxbrab; tnu++) for (var k:Int=0; k<roK; k++) {
                                //Console.OUT.printf("tmu=%d tnu=%d k=%d ind=%d val=%e\n",tmu,tnu,k,ind,temp(ind));
                                val m = norm(tmu)*norm(tnu)*temp(ind++);
                                dk(k) += density(tmu,tnu)*m; // eqn 15b
                                for(aorb in 0..(nOrbital-1)) {
                                    muak(tmu,aorb,k) += mos(aorb,tnu) * m; // eqn 16b the most expensive step!!!
                                }
                                munuk(tmu,tnu,k) = m;
                            }                                
                        else // Be careful... this is tricky ... maxbra are not swap 
                            for (var tnu:Int=nu; tnu<nu+maxbrab; tnu++) for (var tmu:Int=mu; tmu<mu+maxbraa; tmu++) for (var k:Int=0; k<roK; k++) {
                                //Console.OUT.printf("(Swap) tmu=%d tnu=%d k=%d ind=%d val=%e\n",tmu,tnu,k,ind,temp(ind));
                                val m = norm(tmu)*norm(tnu)*temp(ind++);
                                dk(k) += density(tmu,tnu)*m; // eqn 15b
                                for(aorb in 0..(nOrbital-1)) {
                                    muak(tmu,aorb,k) += mos(aorb,tnu) * m; // eqn 16b the most expensive step!!!
                                }
                                munuk(tmu,tnu,k) = m;
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

        for (tmu in 0..(nBasis-1)) for (tnu in 0..(nBasis-1)) {
            var jContrib:Double = 0.0;
            var kContrib:Double = 0.0;
            for (k in 0..(roK-1))  {
                jContrib += munuk(tmu,tnu,k)*dk(k); // eqn 15a
                for (a in 0..(nOrbital-1)) {
                    kContrib += muak(tmu,a,k)*muak(tnu,a,k); // eqn16a 
                }
            }
            jMatrix(tmu,tnu) = jContrib;
            kMatrix(tmu,tnu) = kContrib;
            this(tmu,tnu) = 2.0*jContrib-kContrib;
        }
/*
        Console.OUT.println("diag2E RO");
        Console.OUT.println(diagtwoe);

        Console.OUT.println("J Mat RO");
        Console.OUT.println(jMatrix);

        Console.OUT.println("K Mat RO");
        Console.OUT.println(kMatrix);
*/

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix with RO: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }
 
}

