/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

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

    public def this(mol:Molecule[QMAtom]!,  
                    oneE:OneElectronIntegrals!, 
                    twoE:TwoElectronIntegrals!, gMatType:Int) {
        super(mol, oneE, twoE);

        this.gMatType = gMatType;
    }

    public def scf() : Void {
        // check first if closed shell run?
        val noOfElectrons = molecule.getNumberOfElectrons();
        val noOfOccupancies = noOfElectrons / 2;
        
        if (noOfElectrons%2 != 0) {
           x10.io.Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        } // end if

        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();
        
        energy = 0.0;

        var converged:Boolean = false;
        var oldEnergy:Double = 0.0; 
        var nuclearEnergy:Double = nuclearEnergy();

        x10.io.Console.OUT.println("Nuclear repulsion energy = " + nuclearEnergy + " a.u.");

        var eOne:Double, eTwo:Double;

        Console.OUT.println ("    Initializing matrices ...");
        // init memory for the matrices
        val N = hCore.getRowCount();
        val gMatrix  = new GMatrix(N) as GMatrix!;
        val mos      = new MolecularOrbitals(N) as MolecularOrbitals!;
        val density  = new Density(N) as Density!; // density.make();

        var fock:Fock!  = new Fock(N) as Fock!;

        Console.OUT.println("    Forming initial guess ...");
        // compute initial MOs
        mos.compute(hCore, overlap);

        x10.io.Console.OUT.println("    Starting RHF-SCF ... ");        

        val diis = new DIISFockExtrapolator();

        // start the SCF cycle
        for(var scfIteration:Int=0; scfIteration<maxIteration; scfIteration++) {
            // make or guess density
            density.compute(noOfOccupancies, mos);
            
            // make the G matrix
            gMatrix.compute(twoE, density, gMatType);
           
            val timer = new Timer(2);

            timer.start(0);
            // make fock matrix
            fock.compute(hCore, gMatrix);
            // TODO: DIIS to be switched on / off on request 
            fock = diis.next(fock, overlap, density);
            timer.stop(0);
            Console.OUT.println ("    Time to construct Fock: " + (timer.total(0) as Double) / 1e9 + " seconds");
 
            
            timer.start(1);
            // compute the new MOs
            mos.compute(fock, overlap);
            timer.stop(1);
            Console.OUT.println ("    Time to form MOS: " + (timer.total(1) as Double) / 1e9 + " seconds");
         
            // compute the total energy at this point
            eOne = density.mul(hCore).trace();
            x10.io.Console.OUT.println("    Nuclear electron attraction energy = " + eOne + " a.u.");

            eTwo = density.mul(fock).trace();
            x10.io.Console.OUT.println("    Electron repulsion energy = " + eTwo + " a.u.");
            
            energy = eOne + eTwo + nuclearEnergy;

            x10.io.Console.OUT.println("Cycle #" + scfIteration + " Total energy = " + energy);
                        
            // ckeck for convergence
            if (Math.abs(energy - oldEnergy) < energyTolerance) {
                converged = true;
                break;
            } // end if
            
            oldEnergy = energy;
        } // end of SCF iteration

        if (!converged) 
           x10.io.Console.OUT.println("SCF did not converge in " + maxIteration + " cycles!");
        else 
           x10.io.Console.OUT.println("SCF converged. Final SCF energy = " + energy + " a.u.");
    }
}

