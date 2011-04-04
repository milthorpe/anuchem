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

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class DiatomicForceField implements ForceField {
    val diatomicPotentials : Array[DiatomicPotential](1){rect,zeroBased,rail};

    public def this(diatomicPotentials : Array[DiatomicPotential](1){rect,zeroBased,rail}) {
        this.diatomicPotentials = diatomicPotentials;
    }
    
    public def getPotentialAndForces(atoms: DistArray[Array[MMAtom](1){rect,zeroBased,rail}](1)) : Double {
        var V : Double = 0.0;        
        for (p in 0..(diatomicPotentials.size-1)) {
            val potential = diatomicPotentials(p);
            V += potential.getPotentialAndForces();
        }
        return V;
    }

    public def getAtomMass(symbol : String) : Double {
        if (symbol.equals("H")) {
            return 1.0079;
        } else if (symbol.equals("F")) {
            return 18.9984;
        } else {
            // TODO
            return 10000.0;
        }
    }
}
