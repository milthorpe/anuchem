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

import au.edu.anu.chem.mm.MMAtom;

public abstract class DiatomicPotential {
    public val atom1 : MMAtom;
    public val atom2 : MMAtom;

    public def this(atom1 : MMAtom, atom2 : MMAtom) {
        this.atom1 = atom1;
        this.atom2 = atom2;
    }

    public abstract def getPotentialAndForces() : Double;
}

