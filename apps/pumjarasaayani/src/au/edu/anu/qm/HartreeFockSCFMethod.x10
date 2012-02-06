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

import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;

/**
 * Implementation of Hartree-Fock SCF method
 *
 * @author: V.Ganesh
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
        val noOfElectrons = molecule.getNumberOfElectrons();        
        if (noOfElectrons%2 != 0 || molecule.getMultiplicity()!=1) {
           Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        }

        val nuclearEnergy = getNuclearEnergy();
        Console.OUT.printf("Nuclear repulsion energy = %.6f a.u.\n", nuclearEnergy);

        val noOfBasisFunctions:Long = bfs.getBasisFunctions().size();
        val noOfIntegrals:Long = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;
        Console.OUT.println("\nNumber of 2E integrals: " + noOfIntegrals);

        energy = 0.0;

        var converged:Boolean = false;
        var oldEnergy:Double = 0.0;

        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();

        // init memory for the matrices
        val N = hCore.getRowCount();
        val jd = JobDefaults.getInstance();
        val gMatrixRO  = new GMatrixRO(N, bfs, molecule);
        val gMatrix  = new GMatrix(N, bfs, molecule);

        val mos      = new MolecularOrbitals(N);
        val noOfOccupancies = noOfElectrons / 2;
        val density  = new Density(N, noOfOccupancies); // density.make();

        var fock:Fock  = new Fock(N);

        //Console.OUT.println("    Forming initial guess ...");
   
	if (jd.guess==1) density.applyGuess(bfs.getSAD()); // GUESS = SAD
        else if (jd.guess==0) {         // GUESS = CORE 
            mos.compute(hCore, overlap);
	    density.compute(mos);
        }

        Console.OUT.printf("----------------------------------------------------------\n");      

        val diis = new DIISFockExtrapolator();

        // start the SCF cycle
        for(var scfIteration:Int=0; scfIteration<maxIteration; scfIteration++) {
            // make or guess density
	        if (scfIteration>0) {
                density.compute(mos);
            }
            
            // make the G matrix
            gMatrixRO.compute(density, mos);
            gMatrix.compute(density);

            //Console.OUT.println("G Mat");
            //Console.OUT.println(gMatrix);

            //Console.OUT.println("G Mat RO");
            //Console.OUT.println(gMatrixRO);

            //val timer = new Timer(2);

            //timer.start(0);
            // make fock matrix
            if (jd.roOn>0) fock.compute(hCore, gMatrixRO);
            else fock.compute(hCore, gMatrix);
            // SCF_ALGORITHM = DIIS  
            fock = diis.next(fock, overlap, density);
            //timer.stop(0);
            //Console.OUT.println ("    Time to construct Fock: " + (timer.total(0) as Double) / 1e9 + " seconds");
     
            //timer.start(1);
            // compute the new MOs
            mos.compute(fock, overlap);
            //timer.stop(1);
            //Console.OUT.println ("    Time to form MOS: " + (timer.total(1) as Double) / 1e9 + " seconds");
         
            // compute the total energy at this point
            val eOne = density.mul(hCore).trace();
            // Console.OUT.printf("  Nuclear electron attraction energy = %.6f a.u.\n", eOne);

            val eTwo = density.mul(fock).trace();
            // Console.OUT.printf("  Electron repulsion energy = %.6f a.u.\n", eTwo);
            
            energy = eOne + eTwo + nuclearEnergy;

            Console.OUT.printf("Cycle #%i Total energy = %.6f a.u. (scale factor = %.6f)", scfIteration, energy/roZ,roZ);
            if (scfIteration>0) {
                Console.OUT.printf(" (%.6f)", (energy-oldEnergy)/roZ);
            } else {
                // ignore the first cycle's timings 
                // as far fewer integrals are calculated
                gMatrix.timer.clear(0);
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

        Console.OUT.printf("==========================================================\n");

        Console.OUT.println("GMatrix construction timings:");
        Console.OUT.println("iters     mean   stddev      min      max");
        Console.OUT.printf("%5i %8.4g %8.4g %8.4g %8.4g", 
                            gMatrix.timer.count(0), 
                            (gMatrix.timer.mean(0) as Double) / 1e9,
                            (gMatrix.timer.stdDev(0)) / 1e9,
                            (gMatrix.timer.min(0) as Double) / 1e9,
                            (gMatrix.timer.max(0) as Double) / 1e9);
    }
}

