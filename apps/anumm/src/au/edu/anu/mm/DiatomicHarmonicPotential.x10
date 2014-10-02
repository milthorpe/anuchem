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

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents an Harmonic potential between two atoms
 * of the form V(r) = 1/2 k (r-b)^2 where:
 * r is the bond length;
 * k is the force constant; and
 * b is the equilibrium bond length
 * "In the simplest approximation the potential function is a parabola..." 
 * - Berendsen, "Simulating the Physical World", p.6 (2007)
 */
public class DiatomicHarmonicPotential extends DiatomicPotential {
    /** The force constant in kJ mol^-1 nm^-2. */
    public val forceConstant : Double;

    /** The equilibrium bond length in nm. */
    public val bondLength : Double;

    public def this(atom1 : MMAtom, atom2 : MMAtom, bondLength : Double, forceConstant : Double) {
        super(atom1, atom2);
        this.bondLength = bondLength;
        this.forceConstant = forceConstant;
    }

    public def getPotentialAndForces() : Double {
        val r = atom2.centre - atom1.centre;
        val displacement = r.length() - bondLength;
        Console.OUT.print(r.length() + " ");
        atom1.force = forceConstant * displacement * r.normalize();
        atom2.force = -atom1.force;
        return 0.5 * forceConstant * displacement * displacement;
    }
}

