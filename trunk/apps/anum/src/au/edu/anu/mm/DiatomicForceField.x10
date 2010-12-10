/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class DiatomicForceField implements ForceField {
    val diatomicPotentials : Array[DiatomicPotential](1){rail};

    public def this(diatomicPotentials : Array[DiatomicPotential](1){rail}) {
        this.diatomicPotentials = diatomicPotentials;
    }
    
    public def getPotentialAndForces(atoms: DistArray[Array[MMAtom](1){rail}](1)) : Double {
        var V : Double = 0.0;        
        for([p] in 0..diatomicPotentials.length-1) {
            val potential = diatomicPotentials(p);
            V += potential.getPotentialAndForces();
        }
        return V;
    }
}
