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
package au.edu.anu.mm.uff;

import au.edu.anu.chem.BondType;
import au.edu.anu.mm.LennardJonesParameters;

/**
 * This struct represents the Universal Force Field parameters for given atom type.
 * @property description the atom type e.g. HW-OW-HW for the angle in TIP3P water
 * @property mass the mass of the atom type e.g. 1.00794 for H
 * @property bondLength the equilibrium bond length in Angstroms
 * @property bondAngle the equilibrium bond angle in rad
 * @property vdwParams the Lennard-Jones parameters for non-bond van der Waals interactions
 * @property effectiveCharge the effective charge for use in bond-stretching and angle distortion forces
 * @property electronegativity the GMP electronegativity for bond length calculation
 */
public struct UffParameters(description : String, 
                            mass : Double,
                            bondRadius : Double, 
                            bondAngle : Double,
                            vdwParams : LennardJonesParameters,
                            effectiveCharge : Double,
                            electronegativity : Double) {
    public def this(description : String, mass : Double, bondRadius : Double, bondAngle : Double, vdwParams : LennardJonesParameters, effectiveCharge : Double, electronegativity : Double) {
        property(description, mass, bondRadius, bondAngle, vdwParams, effectiveCharge, electronegativity);
    }
}
