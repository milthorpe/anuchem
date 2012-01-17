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

        val mostMat=mos.getMatrix();
        val denMat=density.getMatrix();
        val nOrbital=desity.getNoOfOccupancies();

        val maxbral = bfs.getShellList().getMaximumAngularMomentum();
        val maxbra = (maxbral+1)*(maxbral+2)/2; 

        val roK = (roN+1)*(roL+1)*(roL+1);
        val dk = new Array[Double](0..(roK-1)); // eqn 15b in RO#7
        val temp = new Array[Double](0..(maxbra*maxbra)*0..(roK-1)); // Result for one batch

        // Infinite memory code
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
               for(var b:Int=0; b<=a; b++) {
                   val bFunc = mol.getAtom(b).getBasisFunctions();
                   // val nbFunc = (b<a) ? bFunc.size() : i+1;
                   // basis functions on b
                   for(var j:Int=0; j<nbFunc; j++) {
                       val jbFunc = bFunc.get(j);

                       // swap A and B if B has higher anglular momentum than A
                       extract info from basisfunctions 
                       call genclass (temp, info from basis function)
                       // transfer infomation from temp to munuk

                       for (mu; ;) for (nu; ;) for (k=0; k<K; k++) 
                           dk(k) += denMat(mu,nu)*munuk(mu,nu,k); // eqn 15b
                       for (mu; ;) for (nu; ;) for (a=0; a<nOrbital; a++) for (k=0; k<K; k++) 
                           muak(mu,a,k) + = mosMat(a,nu) * munuk(mu,nu,k); // eqn 16b the most expensive step!!!

                   }
               }
            }
        }     
        
        for (mu=0; mu<nBasis; mu++) for (nu=0; nu<nBasis; nu++) for (k=0; k<K; k++)  {
            jMat(mu,nu) += munuk(mu.nu,k)*dk(k); // eqn 15a
            for (a=0; a<nOrbital; a++)
                kMat(mu,nu) += muak(mu,a,k)*munuk(nu,a,k); // eqn16a 
        }

        for (mu=0; mu<nBasis; mu++) for (nu=0; nu<nBasis; nu++) 
           gMat(mu,nu) = hCore(mu,nu) + jMat(mu,nu) - 0.5 * kMat(mu,nu); //eqn14

        timer.stop(0);
        Console.OUT.printf("    Time to construct GMatrix: %.3g seconds\n", (timer.last(0) as Double) / 1e9);
    }

 
}

