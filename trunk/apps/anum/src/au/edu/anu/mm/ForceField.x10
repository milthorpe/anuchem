/*
 * This file is part of ANU Molecular Mechanics (ANUMM).
 * ANUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUMM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This interface represents an all-atom force field such as
 * AMBER or the Universal Force Field.
 */
public interface ForceField {
    /**
     * Performs a full calculation of potential and updates forces
     * on each atom in the system.
     * @return the energy of the system // TODO units?
     */
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double;

    /**
     * @return the mass of the given atom type e.g. 1.00794 for H
     */
    public global def getAtomMass(symbol : String) : Double;
}
