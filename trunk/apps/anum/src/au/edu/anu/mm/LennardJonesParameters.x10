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

