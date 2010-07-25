/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a Morse potential between two atoms
 * of the form V(r) = D(1 - e^(-a(r-b)))^2 where:
 * D is the dissociation energy;
 * a is a constant related to the steepness of the curve (here a = 2/b); and
 * b is the equilibrium bond length
 * "The Morse curve is only a convenient analytical expression that has some
 * essential features of a diatomic potential, ... but there is no theoretical
 * justification for this particular form." - Berendsen, "Simulating the Physical World", p.6 (2007)
 * @see P. M. Morse "Diatomic molecules according to the wave mechanics.", Phys. Rev. 34, 57-64 (1929)
 */
public class DiatomicMorsePotential extends DiatomicPotential {
    /** The dissociation energy in kJ mol^-1. */
    public val dissociationEnergy : Double;

    /** The equilibrium bond length in nm. */
    public val bondLength : Double;

    public def this(atom1 : MMAtom, atom2 : MMAtom, bondLength : Double, dissociationEnergy : Double) {
        super(atom1, atom2);
        this.bondLength = bondLength;
        this.dissociationEnergy = dissociationEnergy;
    }

    public def getPotentialAndForces() : Double {
        val r = atom2.centre - atom1.centre;
        val displacement = r.length() - bondLength;
        Console.OUT.println(r.length());
        val a = 2.0 / bondLength;
        val expTerm = Math.exp(-a * displacement);
        // dV/dr = 2D(a*e^(-a(r-b)) (1 - e^(-a(r-b)))
        val expForce = 2.0 * dissociationEnergy * a * expTerm * (1 - expTerm);
        atom1.force = expForce * r.normalize();
        atom2.force = -atom1.force;
        // V(r) = D(1 - e^(-a(r-b)))^2
        return dissociationEnergy * (1.0 - expTerm) * (1.0 - expTerm);
    }
}

