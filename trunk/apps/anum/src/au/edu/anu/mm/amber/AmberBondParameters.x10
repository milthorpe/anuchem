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
package au.edu.anu.mm.amber;

/**
 * This struct represents the AMBER parameters for given bond type.
 * @property description the bond angle type e.g. OW-HW for the bonds in TIP3P water
 * @property forceConstant the bond stretch force constant in kcal/mol*A^2
 * @property bondRadius the equilibrium bond radius in Angstroms
 */
public struct AmberBondParameters(description : String, forceConstant : Double, bondRadius : Double) {
    public def this(description : String, forceConstant : Double, bondRadius : Double) {
        property(description, forceConstant, bondRadius);
    }
}
