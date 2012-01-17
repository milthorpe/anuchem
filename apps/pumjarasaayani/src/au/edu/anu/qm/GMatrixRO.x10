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

    public def this(N:Int, bfs:BasisFunctions, molecule:Molecule[QMAtom]) {
        super(N);
        this.bfs = bfs;
        this.mol = molecule;
        this.nBasis=N;
        val jd = JobDefaults.getInstance();
        this.roN=jd.roN;
        this.roL=jd.roL;
        this.roZ=jd.roZ;  
    }


    public def compute(density:Density, mos:MolecularOrbitals) {
        timer.start(0);

        val N = getRowCount();

        val mosMat = mos.getMatrix();
        val denMat = density.getMatrix();
        val nOrbital = density.getNoOfOccupancies();
        val noOfAtoms = mol.getNumberOfAtoms();

        val roK = (roN+1)*(roL+1)*(roL+1);
        val dk = new Array[Double](0..(roK-1)); // eqn 15b in RO#7
        val K = 1; // TODO

        // Infinite memory code two 3D Arrays
        val munuk = new Array[Double](0..(nBasis-1)*0..(nBasis-1)*0..(roK-1)); // Auxiliary integrals
        val muak = new Array[Double](0..(nBasis-1)*0..(nOrbital-1)*0..(roK-1)); // half-transformed auxiliary integrals eqn 16b in RO#7

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
                    val nbFunc = aFunc.size();
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);

                        // swap A and B if B has higher anglular momentum than A                       

                        // extract info from basisfunctions
                        // Note that iaFunc and jbFunc are ContractedGaussians
                        // val may not work?  must specify the same type as in cpp code?
                        val aang = iaFunc.getTotalAngularMomentum();
                        val aprimitive = iaFunc.getPrimitives();
                        val dConA = aprimitive.size;
                        val apoint = new Rail[Double](3);
                        apoint(0) = aprimitive(0).origin.i;
                        apoint(1) = aprimitive(0).origin.j;
                        apoint(2) = aprimitive(0).origin.k;
                        val conA = new Rail[Double](dConA);
                        val zetaA = new Rail[Double](dConA);
                        for (ai in 0..(dConA-1)) {
                           conA(ai)=aprimitive(ai).coefficient;
                           zetaA(ai)=aprimitive(ai).exponent;
                           conA(ai)*=aprimitive(ai).normalization;
                           // normalization problem to be addressed in integral pack cpp code? -- that will slow things down?
                        }
                        // do the same for b

                        var mu:Int = 1; // TODO
                        var nu:Int = 1; // TODO
                        val bang = 1; // TODO
                        val bpoint = new Rail[Double](3); // TODO
                        val dConB = 1; // TODO
                        val conB = new Rail[Double](dConB); // TODO
                        val zetaB = new Rail[Double](dConB); // TODO
                        val L = 1; // TODO

                        val maxbraa = (aang+1)*(aang+2)/2; 
                        val maxbrab = (bang+1)*(bang+2)/2; 
                        val temp = new Array[Double](0..(maxbraa*maxbrab)*0..(roK-1)); // Result for one batch

                        // call genclass (temp, info from basis function)            
                        // roN, roL should be there during initialization...
                        val aux = new Integral_Pack(N);
                        aux.genClass(aang, bang, apoint, bpoint, zetaA, zetaB, conA, conB, dConA, dConB, L/*, temp*/);      

                        // transfer infomation from temp to munuk (Swap A and B again if necessary)

                        for (tmu in mu..(mu+maxbraa-1)) for (tnu in 0..(nu+maxbrab-1)) for (k in 0..(K-1)) 
                           dk(k) += denMat(tmu,tnu)*munuk(tmu,tnu,k); // eqn 15b
                        for (tmu in mu..(mu+maxbraa-1)) for (tnu in 0..(nu+maxbrab-1)) for (aa in 0..(nOrbital-1)) for (k in 0..(K-1)) 
                           muak(tmu,aa,k) += mosMat(aa,tnu) * munuk(tmu,tnu,k); // eqn 16b the most expensive step!!!

                        mu+=maxbraa; nu+=maxbrab;

                    }
                }
            }
        }     

        val jMatrix = new Matrix(N);
        val kMatrix = new Matrix(N);
        val jMat = jMatrix.getMatrix();
        val kMat = kMatrix.getMatrix();
        val gMat = getMatrix();
        
        for (mu in 0..(nBasis-1)) for (nu in 0..(nBasis-1)) for (k in 0..(K-1))  {
            jMat(mu,nu) += munuk(mu,nu,k)*dk(k); // eqn 15a
            for (a in 0..(nOrbital-1))
                kMat(mu,nu) += muak(mu,a,k)*munuk(nu,a,k); // eqn16a 
        }

        for (mu in 0..(nBasis-1)) for (nu in 0..(nBasis-1))
           gMat(mu,nu) = jMat(mu,nu) - 0.5 * kMat(mu,nu); // eqn14

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }
 
}

