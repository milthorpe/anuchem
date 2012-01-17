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
 * G matrix in HF calculation
 * Integral screening by cutoff based on Häser & Ahlrichs (1989)
 * @see Häser, M. and Ahlrichs, R. (1989). "Improvements on the Direct SCF
 *   method". J. Comp. Chem. 10 (1) pp.104-111.
 *
 * @author: V.Ganesh
 */
public class GMatrixRO extends Matrix {
    public static DEFAULT_GMATTYPE=4;

    public val timer = new StatisticalTimer(1);
    public static TIMER_IDX_TOTAL = 0;

    private val gMatType : Int;
    private val computeInst : DistArray[ComputePlace](1);

    private val bfs : BasisFunctions;
    private val mol : Molecule[QMAtom];

    static THRESH:Double = 1.0e-8;
    private var thresh2:Double = 0.0;

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

        this.gMatType = jd.gMatrixParallelScheme;
        // TO DO   

        val dCut = new Matrix(N);
        val qCut = new Matrix(N);
        // Schwarz cut-off: Häser & Ahlrichs eqn 12
        val shellList = bfs.getShellList();
        val maxam = shellList.getMaximumAngularMomentum();
        val twoE = new TwoElectronIntegrals(maxam, bfs.getNormalizationFactors(), THRESH);


        val fakeDensity = new Density(N, 2);
        val fakeD = fakeDensity.getMatrix();
        fakeD.fill(1.0);

        val noOfAtoms = mol.getNumberOfAtoms();
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
                    val nbFunc = (b<a) ? bFunc.size() : i+1;
                    // basis functions on b
                    for(var j:Int=0; j<nbFunc; j++) {
                        val jbFunc = bFunc.get(j);
                        twoE.compute2EAndRecord(iaFunc, jbFunc, iaFunc, jbFunc, shellList, qCut, dCut, fakeDensity);
                    }
                }
            }
        }

        val est = qCut.getMatrix();

        var maxEst:Double = 0.0;
        for(var a:Int=0; a<N; a++) for(var b:Int=0; b<N; b++) {
            if (a==b) est(a,b)*=.5; else est(a,b)*=.25;
            est(a,b)=Math.sqrt(Math.abs(est(a,b)));
            if (est(a,b)>maxEst) maxEst=est(a,b);
            //Console.OUT.printf("%d %d %e \n",a,b,EST(a,b));
        }

        Console.OUT.printf("\tmaxEst %.4g\n", maxEst);
        val maxEstVal = maxEst;

        computeInst = DistArray.make[ComputePlace](Dist.makeUnique(), (Point) => new ComputePlace(N, molecule, bfs, qCut, dCut, maxEstVal));
    }


    public def compute(density:Density, mos:MolecularOrbitals) {
        timer.start(0);

        val mostMat=mos.getMatrix();
        val denMat=density.getMatrix();
        val nOrbital=mos.something;

        // Infinite memory code
        val maxbra = 6; // dfunctions
        val roK = (roN+1)*(roL+1)*(roL+1);
        val dk = new Array[Double](0..(roK-1)); // eqn 15b in RO#7
        val munuk = new Array[Double](0..(nBasis-1)*0..(nBasis-1)*0..(roK-1)); // Auxiliary integrals
        val muak = new Array[Double](0..(nBasis-1)*0..(nOrbital-1)*0..(roK-1)); // half-transformed auxiliary integrals eqn 16b in RO#7
        val temp = new Array[Double](0..(maxbra*maxbra)*0..(roK-1)); // Result for one batch

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

