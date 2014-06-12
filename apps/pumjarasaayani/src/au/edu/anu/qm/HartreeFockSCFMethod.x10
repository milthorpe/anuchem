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

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;

import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;
import au.edu.anu.util.StatisticalTimer;

import x10.matrix.DenseMatrix;

/**
 * Implementation of Hartree-Fock SCF method
 * @author: V. Ganesh, J. Milthorpe, T. Limpanuparb 
 */
public class HartreeFockSCFMethod extends SCFMethod { 
    val roZ:Double;

    public def this(mol:Molecule[QMAtom],  
                    oneE:OneElectronIntegrals, 
                    bfs:BasisFunctions) {
        super(mol, oneE, bfs);
        val jd = JobDefaults.getInstance();
        this.roZ=jd.roZ;
    }

    public def scf() : void {
        val TIMER_TOTAL:Int = 0n; // Can't use static
        val TIMER_2:Int = 1n;
        val TIMER_S:Int = 2n;
        val TIMER_F:Int = 3n;
        val TIMER_M:Int = 4n;
        val TIMER_GUESS:Int = 5n;
        val timer = new StatisticalTimer(6);

        val noOfElectrons = molecule.getNumberOfElectrons();        
        if (noOfElectrons%2n != 0n || molecule.getMultiplicity()!=1n) {
           Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        }
        val noOfOccupancies = noOfElectrons / 2;

        val nuclearEnergy = getNuclearEnergy();
        Console.OUT.printf("Nuclear repulsion energy = %.10f a.u.\n", nuclearEnergy/roZ);

        /*val noOfBasisFunctions:Long = bfs.getBasisFunctions().size();
        val noOfIntegrals:Long = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;
        Console.OUT.println("\nNumber of 2E integrals: " + noOfIntegrals);*/

        energy = 0.0;

        var converged:Boolean = false;
        var oldEnergy:Double = 0.0;
        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();
        @Ifdef("__DEBUG__") {
            // See3.4.5 Orthogonalization of the Basis pp.142-145 SO Book
            Console.OUT.printf("Checking the overlap matrix...\n"); 
            val diag = new GMLDiagonalizer();
            diag.diagonalize(overlap);
            val eigval =  diag.getEigenValues().d;;
            Console.OUT.printf("Three smallest eigenvalues %.3e %.3e %.3e\n",eigval(0),eigval(1),eigval(2));        
            if (eigval(0)<1e-6) Console.OUT.printf("Linear dependence detected!!!\n");            
        }
        // init memory for the matrices
        val N = hCore.N;
        val jd = JobDefaults.getInstance();

        val gMatrix:GMatrix{self.N==N};
        val roFockMethod:ROFockMethod; // change version here 1st out of 3 places
        val gMatrixRo:DenseMatrix(N,N);
        val roThresh=jd.roThresh;
        val thresh=jd.thresh;
        if (jd.roOn > 0n && maxIteration > 0n) {
            // change version here 2nd out of 3 places
            roFockMethod = new ROFockMethod(N, bfs, molecule, noOfOccupancies, 0., roZ*roThresh);
            gMatrixRo = new DenseMatrix(N,N);
        } else {
            roFockMethod = null;
            gMatrixRo = null;
        }
        if ((jd.roOn==0n || jd.compareRo==true) && maxIteration>0n) {
            gMatrix = new GMatrix(N, bfs, molecule,0.,roZ*thresh);
        } else {
            gMatrix = null;
        }

        val mos = new MolecularOrbitals(N);
        val density = new Density(N, noOfOccupancies);

        var fock:Fock{self.M==N,self.N==N} = new Fock(N);

        timer.start(TIMER_GUESS);
        Console.OUT.printf("Making guess orbitals ");
        if (jd.guess.equals(JobDefaults.GUESS_SAD)) {
            Console.OUT.printf(" (density = SAD, MOs = core)...\n");
            density.applyGuess(bfs.getSAD());
            mos.compute(hCore, overlap); // Cheat
        } else if (jd.guess.equals(JobDefaults.GUESS_CORE)) {
            Console.OUT.printf(" (MOs & density = core)...\n");
            mos.compute(hCore, overlap);
            density.compute(mos);
        } else {
            Console.OUT.printf("(no guess)...\n");
        }
        timer.stop(TIMER_GUESS);
        Console.OUT.println ("    Time to construct guess density & MOs: " + (timer.total(TIMER_GUESS) as Double) / 1e9 + " seconds");
        
        if (maxIteration > 0n) {
            Console.OUT.printf("Starting SCF procedure...\n");      

            val diis = new DIISFockExtrapolator(N);

            // start the SCF cycle
            for(var scfIteration:Int=0n; scfIteration<maxIteration; scfIteration++) {
                // make or guess density
	            if (scfIteration>0) {
                    density.compute(mos);
                }
                
                // make the G matrix
                if (jd.roOn>0) {
                    roFockMethod.compute(density, mos, gMatrixRo);
                    if (jd.compareRo==true) {
                        gMatrix.compute(density);
                        Console.OUT.println("G Mat");
                        Console.OUT.println(gMatrix);

                        Console.OUT.println("G Mat RO");
                        Console.OUT.println(gMatrixRo);
                    }
                } else {
                    gMatrix.compute(density);
                }            

                timer.start(TIMER_F);
                // make fock matrix
                if (jd.roOn>0) fock.compute(hCore, gMatrixRo);
                else fock.compute(hCore, gMatrix);
                // SCF_ALGORITHM = DIIS  
                fock = diis.next(fock, overlap, density);
                timer.stop(TIMER_F);
                Console.OUT.println ("    Time to construct Fock: " + (timer.total(TIMER_F) as Double) / 1e9 + " seconds");

                timer.start(TIMER_M);
                // compute the new MOs
                mos.compute(fock, overlap);
                timer.stop(TIMER_M);
                Console.OUT.println ("    Time to form MOS: " + (timer.total(TIMER_M) as Double) / 1e9 + " seconds");
             
                // compute the total energy 
                val eOne = .5*density.clone().mult(density, hCore).trace();
                val eTwo = .5*density.clone().mult(density, fock).trace();           
                energy = eOne + eTwo + nuclearEnergy;
                @Ifdef("__DEBUG__") { Console.OUT.printf("eOne = %.10f eTwo= %.10f\n",eOne/roZ,eTwo/roZ);}
                Console.OUT.printf("Cycle #%i Total energy = %.10f\n", scfIteration, energy/roZ);
                if (scfIteration>0) {
                    Console.OUT.printf(" (%.6f)", (energy-oldEnergy)/roZ);
                } else {
                    // ignore the first cycle's timings 
                    // as far fewer integrals are calculated
                    if (jd.roOn==0n || jd.compareRo==true) {
                        gMatrix.timer.clear(0);
                    }
                    if (jd.roOn>0n) {
                        roFockMethod.timer.clear(0);
                    }
                }
                Console.OUT.printf("\n----------------------------------------------------------\n");
                // check for convergence
                if (scfIteration > 0 && diis.isConverged()) {
                    converged = true;
                    break;
                } // end if
                
                oldEnergy = energy;
            } // end of SCF iteration

            if (!converged) 
               Console.OUT.println("SCF did not converge in " + maxIteration + " cycles!");
            else 
               Console.OUT.printf("SCF converged. Final SCF energy = %.6f a.u.\n", energy/roZ,roZ);

            if ((jd.roOn == 0n || jd.compareRo)) {
                Console.OUT.println("GMatrix construction timings:");
                gMatrix.timer.printSeconds();
            }

            if (jd.roOn>0n) {
                Console.OUT.println("ROFockMethod construction timings:");
                roFockMethod.timer.printSeconds();
            }
        }

        Console.OUT.printf("==========================================================\n");

        var ecount:Double=0.; var maxDen:Double=0.; var minDen:Double=1e100;
        for (var i:Int=0n; i<N; i++) for (var j:Int=i; j<N; j++) {
             val contrib=density(i,j)*overlap(i,j);
             ecount+=contrib;
             if (i!=j) ecount+=contrib;
             if (density(i,j)>100.) Console.OUT.printf("density(%d,%d)=%e, s=%e\n",i,j,density(i,j),overlap(i,j));
             if (density(i,j)>maxDen) maxDen=density(i,j);
             else if (density(i,j)<minDen) minDen=density(i,j);
        }
        Console.OUT.printf("ecount=%e maxDen=%e minDen=%e\n",ecount,maxDen,minDen);
    
        // long range energy
        //Console.OUT.println("before RO heapSize = " + System.heapSize());
        Console.OUT.printf("Long-range - RO\n");
        val roFockMethodL = new ROFockMethod(N, bfs, molecule, noOfOccupancies, jd.roZ*jd.omega, roZ*jd.roThresh); // RO Thesis Eq (2.22) // change version here 3rd out of 3 places
        val gMatrixRoL = new DenseMatrix(N, N);
        roFockMethodL.compute(density, mos, gMatrixRoL);

        if (jd.compareRo) {
            Console.OUT.printf("Long-range - Conventional\n");
            val gMatrixL = new GMatrix(N, bfs, molecule,jd.roZ*jd.omega,roZ*jd.thresh); // RO Thesis Eq (2.22)
            gMatrixL.compute(density);   
            /*Console.OUT.printf("j=%.5e %.5e k=%.5e %.5e\n",gMatrixL.jMatrix(0,0),gMatrixRoL.jMatrix(0,0),gMatrixL.kMatrix(0,0),gMatrixRoL.kMatrix(0,0));
            var jrms:Double=0.,krms:Double=0.,jmax:Double=0.,kmax:Double=0.;
            for ([x,y] in 0..(N-1)*0..(N-1)) {
                val dj=Math.abs(.5*gMatrixL.jMatrix(x,y)-gMatrixRoL.jMatrix(x,y));
                val dk=Math.abs(.25*gMatrixL.kMatrix(x,y)-2.*gMatrixRoL.kMatrix(x,y));
                jrms=dj*dj;
                krms=dk*dk;
                if (dj>jmax) jmax=dj;
                if (dk>kmax) kmax=dk;
            }
            val NN=N*N;
            jrms=Math.pow(jrms/NN,.5);
            krms=Math.pow(krms/NN,.5);
            Console.OUT.printf("j=%.5e (%.5e) k=%.5e (%.5e)\n",jrms,jmax,krms,kmax);
            Console.OUT.printf("log: j=%.5e (%.5e) k=%.5e (%.5e)\n",-Math.log10(jrms),-Math.log10(jmax),-Math.log10(krms),-Math.log10(kmax));    */   
        }
        //Console.OUT.println("after conventional = " + System.heapSize());
        //fock.compute(hCore, gMatrixRo);
        //val eOne = density.clone().mult(density, hCore).trace();
        //val eTwo = density.clone().mult(density, fock).trace();
        //energy = eOne + eTwo + nuclearEnergy;
        //Console.OUT.printf("Cycle ** Total energy = %.6f a.u. (scale factor = %.6f)",  energy/roZ,roZ);
    }

}

