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

public abstract class DiatomicPotential {
    public val atom1 : MMAtom;
    public val atom2 : MMAtom;

    public def this(atom1 : MMAtom, atom2 : MMAtom) {
        this.atom1 = atom1;
        this.atom2 = atom2;
    }

    public abstract def getPotentialAndForces() : Double;
}

