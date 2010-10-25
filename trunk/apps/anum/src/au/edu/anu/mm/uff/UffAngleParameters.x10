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

/**
 * This struct represents the AMBER parameters for given bond angle type.
 * @property description the bond angle type e.g. HW-OW-HW for the angle in TIP3P water
 * @property forceConstant the angle distortion force constant in kcal/mol*rad^2
 * @property bondAngle the equilibrium bond angle in rad
 */
public struct AmberAngleParameters(description : String, forceConstant : Double, bondAngle : Double) {
    public def this(description : String, forceConstant : Double, bondAngle : Double) {
        property(description, forceConstant, bondAngle);
    }
}
