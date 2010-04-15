/**
 * SCFMethod.x10
 *
 * Stub for SCF method
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

public abstract class SCFMethod { 
    protected global val molecule:Molecule[QMAtom]{self.at(this)};
    protected global val oneE:OneElectronIntegrals{self.at(this)};
    protected global val twoE:TwoElectronIntegrals{self.at(this)};

    protected var maxIteration:Int;
    protected var energyTolerance:Double;

    protected var energy:Double;

    public def this(mol:Molecule[QMAtom]!,  
                    oneE:OneElectronIntegrals!, 
                    twoE:TwoElectronIntegrals!) { 
        this.molecule = mol;
        this.oneE     = oneE;
        this.twoE     = twoE;

        val jd = JobDefaults.getInstance();

        maxIteration    = at(jd) { jd.getMaxIterations() };
        energyTolerance = at(jd) { jd.getEnergyTolerance() }; 

        energy = 0.0;
    } 

    public abstract def scf() : Void;

    public def getMaxIteration() = maxIteration;
    public def setMaxIteration(mxIter:Int) : Void {
        maxIteration = mxIter;
    }

    public def getEnergyTolerance() = energyTolerance;
    public def setEnergyTolerance(eneTol:Double) : Void {
        energyTolerance = eneTol;
    }

    public def getEnergy() = energy;

    public def nuclearEnergy() : Double {
        var eNuke:Double = 0.0;
        var i:Int, j:Int;
        val noOfAtoms = molecule.getNumberOfAtoms();
        
        var atomI:QMAtom{self.at(this)}, atomJ:QMAtom{self.at(this)};
        
        // read in the atomic numbers
        atomicNumbers:Array[Int](1)! = new Array[Int]([0..noOfAtoms]);
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
                         / atomI.centre.distance(atomJ.centre);
            } // end for
        } // end for
        
        return eNuke;
    }
}

