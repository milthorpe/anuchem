/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class PenningForceField implements ForceField {

    public val magneticField:Vector3d;
    // static electricField = ??

    public def this(B:Vector3d) {
        this.magneticField = B;
    }
    
    public def getPotentialAndForces(atoms: DistArray[Rail[MMAtom]](1)) : Double {
        var V : Double = 0.0;
        finish ateach(place in atoms) {
            val atomsHere = atoms(place);
            for ([p] in atomsHere) {
                val atom = atomsHere(p);
                //Console.OUT.println("atom.velocity = " + atom.velocity);
                val F = atom.charge * atom.velocity.cross(magneticField);
                //Console.OUT.println("F = " + F);
                atom.force = F;
            }
        }
        return V;
    }

    public def getAtomMass(symbol : String) : Double {
        if (symbol.equals("H")) {
            return 1.0079;
        } else if (symbol.equals("F")) {
            return 18.9984;
        } else {
            throw new IllegalArgumentException("no atom mass found for symbol " + symbol);
        }
    }
}
