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
        val noOfOccupancies = noOfElectrons / 2;

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
        val N = hCore.N;
        val jd = JobDefaults.getInstance();
        val gMatrixRo:GMatrixROmem{self.N==N};
        if (jd.roOn>0) {
            gMatrixRo = new GMatrixROmem(N, bfs, molecule, noOfOccupancies);
        } else {
            gMatrixRo = null;
        }
        val gMatrix = new GMatrix(N, bfs, molecule);

        val mos = new MolecularOrbitals(N);
        val density = new Density(N, noOfOccupancies); // density.make();

        var fock:Fock{self.M==N,self.N==N} = new Fock(N);

        if (jd.guess.equals(JobDefaults.GUESS_SAD)) {
            density.applyGuess(bfs.getSAD());
        } else { // if (jd.guess.equals(JobDefaults.GUESS_CORE)) {
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
            if (jd.roOn>0) {
                gMatrixRo.compute(density, mos);

                if (jd.compareRo) {
                    gMatrix.compute(density);
                    Console.OUT.println("G Mat");
                    Console.OUT.println(gMatrix);

                    Console.OUT.println("G Mat RO");
                    Console.OUT.println(gMatrixRo);
                }
            } else {
                gMatrix.compute(density);
            }

            //val timer = new Timer(2);

            //timer.start(0);
            // make fock matrix
            if (jd.roOn>0) fock.compute(hCore, gMatrixRo);
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
            val eOne = density.clone().mult(density, hCore).trace();
            // Console.OUT.printf("  Nuclear electron attraction energy = %.6f a.u.\n", eOne);

            val eTwo = density.clone().mult(density, fock).trace();
            // Console.OUT.printf("  Electron repulsion energy = %.6f a.u.\n", eTwo);
            
            energy = eOne + eTwo + nuclearEnergy;

            Console.OUT.printf("Cycle #%i Total energy = %.6f a.u. (scale factor = %.6f)", scfIteration, energy/roZ,roZ);
            if (scfIteration>0) {
                Console.OUT.printf(" (%.6f)", (energy-oldEnergy)/roZ);
            } else {
                // ignore the first cycle's timings 
                // as far fewer integrals are calculated
                gMatrix.timer.clear(0);
                if (jd.roOn>0) {
                    gMatrixRo.timer.clear(0);
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

        Console.OUT.printf("==========================================================\n");

        if (jd.roOn == 0 || jd.compareRo) {
            Console.OUT.println("GMatrix construction timings:");
            gMatrix.timer.printSeconds();
        }

        if (jd.roOn>0) {
            Console.OUT.println("GMatrixROmem construction timings:");
            gMatrixRo.timer.printSeconds();
        }
    }
}

