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

/**
 * Represents the parameters for a Lennard-Jones potential, e.g. to 
 * model the van der Waals interaction between non-bonded particles.
 * @property scale the scale parameter i.e. the exponent of the repulsion term
 * @property vdwBondRadius the van der Waals bond radius or equilibrium separation distance
 * @property vdwWellDepth the van der Waals well depth or minimum energy in hartrees in Angstroms
 */
public struct LennardJonesParameters(description : String, 
                                   scale : Double,
                                   equlibrium : Double,
                                   wellDepth : Double) {
    public def this(description : String, scale : Double, bondRadius : Double, wellDepth : Double) {
        property(description, scale, bondRadius, wellDepth);
    }
}

