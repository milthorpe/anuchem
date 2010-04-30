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
 * This struct represents the AMBER parameters for given bond angle type.
 * @property description the bond angle type e.g. HW-OW-HW for the angle in TIP3P water
 * @property forceConstant the angle distortion force constant in kcal/mol*rad^2
 * @property bondAngle the equilibrium bond angle in rad
 */
public struct AmberAngleParameters(description : String, forceConstant : Double, bondAngle : Double) {
    public this(description : String, forceConstant : Double, bondAngle : Double) {
        property(description, forceConstant, bondAngle);
    }
}
