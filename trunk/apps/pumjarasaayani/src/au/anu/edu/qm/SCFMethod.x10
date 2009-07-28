/**
 * SCFMethod.x10
 *
 * Stub for SCF method
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

public abstract class SCFMethod { 
    protected var molecule:Molecule;
    protected var oneE:OneElectronIntegrals;
    protected var twoE:TwoElectronIntegrals;

    protected var maxIteration:Int;
    protected var energyTolerance:Double;

    public def this(mol:Molecule, oneE:OneElectronIntegrals, twoE:TwoElectronIntegrals) { 
        this.molecule = mol;
        this.oneE     = oneE;
        this.twoE     = twoE;

        maxIteration    = 100;
        energyTolerance = 1.0e-5; 
    } 

    public abstract def scf() : void;

    public def getMaxIteration() : Int = maxIteration;
    public def setMaxIteration(mxIter:Int) : void {
        maxIteration = mxIter;
    }

    public def getEnergyTolerance() : Double = energyTolerance;
    public def setEnergyTolerance(eneTol:Double) : void {
        energyTolerance = eneTol;
    }

    public def nuclearEnergy() : Double {
        var eNuke:Double = 0.0;
        var i:Int, j:Int;
        val noOfAtoms = molecule.getNumberOfAtoms();
        
        var atomI:Atom, atomJ:Atom;
        
        // read in the atomic numbers
        atomicNumbers:Array[Int]{rank==1} = Array.make[Int]([0..noOfAtoms]);
        val ai = AtomInfo.getInstance();
        
        for(i=0; i<noOfAtoms; i++) {
            atomicNumbers(i) = ai.getAtomicNumber(molecule.getAtom(i));
        } // end for
        
        // and compute nuclear energy
        for(i=0; i<noOfAtoms; i++) {
            atomI = molecule.getAtom(i);
            for(j=0; j<i; j++) {
                atomJ = molecule.getAtom(j);
                
                eNuke += atomicNumbers(i) * atomicNumbers(j) 
                         / atomI.distanceFrom(atomJ);
            } // end for
        } // end for
        
        return eNuke;
    }
}

