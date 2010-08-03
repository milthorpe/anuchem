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
