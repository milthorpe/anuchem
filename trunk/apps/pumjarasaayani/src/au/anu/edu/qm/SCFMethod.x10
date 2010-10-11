/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import au.edu.anu.chem.AtomInfo;
import au.edu.anu.chem.Molecule;

/**
 * SCFMethod.x10
 *
 * Stub for SCF method
 *
 * @author: V.Ganesh
 */
public abstract class SCFMethod { 
    protected val molecule:Molecule[QMAtom];
    protected val oneE:OneElectronIntegrals;
    protected val bfs:BasisFunctions;

    protected var maxIteration:Int;
    protected var energyTolerance:Double;

    protected var energy:Double;

    public def this(mol:Molecule[QMAtom],  
                    oneE:OneElectronIntegrals,
                    bfs:BasisFunctions) { 
        this.molecule = mol;
        this.oneE = oneE;
        this.bfs = bfs;

        val jd = JobDefaults.getInstance();

        maxIteration    = jd.getMaxIterations();
        energyTolerance = jd.getEnergyTolerance(); 

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

    public def getNuclearEnergy() : Double {
        var eNuke:Double = 0.0;
        var i:Int, j:Int;
        val noOfAtoms = molecule.getNumberOfAtoms();
        
        var atomI:QMAtom, atomJ:QMAtom;
        
        // read in the atomic numbers
        val atomicNumbers = new Array[Int](noOfAtoms);
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

