/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2011.
 */
package au.edu.anu.mm;

import x10.regionarray.DistArray;
import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class DiatomicForceField implements ForceField {
    public static SPECIES_H = 0n;
    public static SPECIES_F = 1n;

    val diatomicPotentials : Rail[DiatomicPotential];

    public def this(diatomicPotentials : Rail[DiatomicPotential]) {
        this.diatomicPotentials = diatomicPotentials;
    }
    
    public def computePotentialAndForces(atoms:DistArray[Rail[MMAtom]](1)):Double {
        var V : Double = 0.0;        
        for (p in 0..(diatomicPotentials.size-1)) {
            val potential = diatomicPotentials(p);
            V += potential.getPotentialAndForces();
        }
        return V;
    }

    public def getAtomMass(species:Int) : Double {
        switch(species) {
            case SPECIES_H :
                return 1.0079;
            case SPECIES_F :
                return 18.9984;
            default :
                throw new IllegalArgumentException("no atom mass found for species " + species);
        }
    }
}
