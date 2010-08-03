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
