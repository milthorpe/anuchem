/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import au.edu.anu.chem.Molecule;
import au.edu.anu.util.Timer;

/**
 * HartreeFockSCFMethod.x10
 *
 * Implementation of Hartree-Fock SCF method
 *
 * @author: V.Ganesh
 */
public class HartreeFockSCFMethod extends SCFMethod { 
    val gMatType:Int;

    public def this(mol:Molecule[QMAtom],  
                    oneE:OneElectronIntegrals, 
                    bfs:BasisFunctions, gMatType:Int) {
        super(mol, oneE, bfs);

        this.gMatType = gMatType;
    }

    public def scf() : void {
        val noOfBasisFunctions = bfs.getBasisFunctions().size();
        val noOfIntegrals:Long = noOfBasisFunctions * (noOfBasisFunctions + 1)
                          * (noOfBasisFunctions * noOfBasisFunctions
                             + noOfBasisFunctions + 2) / 8;
        Console.OUT.println("\nNumber of 2E integrals: " + noOfIntegrals);
    
        // check first if closed shell run?
        val noOfElectrons = molecule.getNumberOfElectrons();
        val noOfOccupancies = noOfElectrons / 2;
        
        if (noOfElectrons%2 != 0) {
           Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        } // end if

        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();
        
        energy = 0.0;

        var converged:Boolean = false;
        var oldEnergy:Double = 0.0;

        val nuclearEnergy = getNuclearEnergy();
        Console.OUT.printf("Nuclear repulsion energy = %.6f a.u.\n", nuclearEnergy);

        //Console.OUT.println ("    Initializing matrices ...");
        // init memory for the matrices
        val N = hCore.getRowCount();
        val gMatrix  = new GMatrix(N, bfs, molecule, gMatType);
        val mos      = new MolecularOrbitals(N);
        val density  = new Density(N, noOfOccupancies); // density.make();

        var fock:Fock  = new Fock(N);

        //Console.OUT.println("    Forming initial guess ...");
        // compute initial MOs
        mos.compute(hCore, overlap);
	    density.compute(mos);
	    density.applyGuess(bfs.getSAD());

        //Console.OUT.println("    Starting RHF-SCF ... ");        

        var diis:DIISFockExtrapolator = null;

        // start the SCF cycle
        for(var scfIteration:Int=0; scfIteration<maxIteration; scfIteration++) {
            // make or guess density
	    if (scfIteration>0) //
            density.compute(mos);
            
            // make the G matrix
            gMatrix.compute(density);
           
            val timer = new Timer(2);

            timer.start(0);
            // make fock matrix
            fock.compute(hCore, gMatrix);
            // TODO: DIIS to be switched on / off on request 
            if (scfIteration % 50==0) diis = new DIISFockExtrapolator();
            fock = diis.next(fock, overlap, density);
            timer.stop(0);
            //Console.OUT.println ("    Time to construct Fock: " + (timer.total(0) as Double) / 1e9 + " seconds");
 
            
            timer.start(1);
            // compute the new MOs
            mos.compute(fock, overlap);
            timer.stop(1);
            //Console.OUT.println ("    Time to form MOS: " + (timer.total(1) as Double) / 1e9 + " seconds");
         
            // compute the total energy at this point
            val eOne = density.mul(hCore).trace();
            Console.OUT.printf("    Nuclear electron attraction energy = %.6f a.u.\n", eOne);

            val eTwo = density.mul(fock).trace();
            Console.OUT.printf("    Electron repulsion energy = %.6f a.u.\n", eTwo);
            
            energy = eOne + eTwo + nuclearEnergy;

            Console.OUT.printf("Cycle #" + scfIteration + " Total energy = %.6f a.u.", energy);
            if (scfIteration>0) Console.OUT.printf(" (%.6f)",energy-oldEnergy);
            Console.OUT.printf("\n");
            // ckeck for convergence
            if (scfIteration>0/*Math.abs(energy - oldEnergy) < energyTolerance &&*/&& diis.isConverged()) {
                converged = true;
                break;
            } // end if
            
            oldEnergy = energy;
        } // end of SCF iteration

        if (!converged) 
           Console.OUT.println("SCF did not converge in " + maxIteration + " cycles!");
        else 
           Console.OUT.printf("SCF converged. Final SCF energy = %.6f a.u.\n", energy);
    }
}

