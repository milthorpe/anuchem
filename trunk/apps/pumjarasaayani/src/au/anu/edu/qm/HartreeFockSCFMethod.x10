/**
 * HartreeFockSCFMethod.x10
 *
 * Implementation of Hartree-Fock SCF method 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import au.edu.anu.chem.Molecule;

public class HartreeFockSCFMethod extends SCFMethod { 
    public def this(mol:Molecule[QMAtom]!,  
                    oneE:OneElectronIntegrals!, 
                    twoE:TwoElectronIntegrals!) {
        super(mol, oneE, twoE);
    }

    public def scf() : void {
        x10.io.Console.OUT.println("Starting RHF-SCF ... ");        
        
        // check first if closed shell run?
        val noOfElectrons = molecule.getNumberOfElectrons();
        val noOfOccupancies = noOfElectrons / 2;
        
        if (noOfElectrons%2 != 0) {
           x10.io.Console.OUT.println("Currently supports only closed shell calculations!");
           return;
        } // end if

        val hCore   = oneE.getHCore();
        val overlap = oneE.getOverlap();
        
        var converged:Boolean = false;
        var energy:Double = 0.0;
        var oldEnergy:Double = 0.0; 
        var nuclearEnergy:Double = nuclearEnergy();

        x10.io.Console.OUT.println("Nuclear repulsion energy = " + nuclearEnergy + " a.u.");

        var eOne:Double, eTwo:Double;

        // init memory for the matrices
        val N = hCore.getRowCount();
        val gMatrix:GMatrix{self.at(this)}       = GMatrix.make(new GMatrix(), N) as GMatrix{self.at(this)};
        val mos:MolecularOrbitals{self.at(this)} = MolecularOrbitals.make(new MolecularOrbitals(), N) as MolecularOrbitals{self.at(this)};
        val density:Density{self.at(this)}       = Density.make(new Density(), N) as Density{self.at(this)};
        val fock:Fock{self.at(this)}             = Fock.make(new Fock(), N) as Fock{self.at(this)};

        // compute initial MOs
        mos.compute(hCore, overlap);

        // start the SCF cycle
        for(var scfIteration:Int=0; scfIteration<maxIteration; scfIteration++) {
            // make or guess density
            density.compute(noOfOccupancies, mos);
            
            // make the G matrix
            gMatrix.compute(twoE, density);
           
            val timer = new Timer(2);

            timer.start(0);
            // make fock matrix
            fock.compute(hCore, gMatrix);
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

