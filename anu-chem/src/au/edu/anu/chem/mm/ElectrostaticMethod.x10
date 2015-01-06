/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright IBM Corporation 2014.
 */

package au.edu.anu.chem.mm;

/**
 * Provides a method for calculating electrostatic interactions between
 * pairs of charged particles in a molecular mechanics simulation.
 */
public interface ElectrostaticMethod {
    /**
     * Compute potential and forces for all particles.
     * @return total system potential
     */
    public def computePotentialAndForces():Double;
}

