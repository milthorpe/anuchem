/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Josh Milthorpe 2010.
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.regionarray.DistArray;
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
    public def getPotentialAndForces(particleDataPlh:PlaceLocalHandle[ParticleData]):Double;

    /**
     * @return the mass of the given atom type e.g. 1.00794 for H=1
     */
    public def getAtomMass(species:Int) : Double;

    /**
     * @return the list of all species supported by this force field
     */
    public def getAtomTypes():Rail[AtomType];

    /**
     * @return the bond type index within this force field for the given species
     */
    public def getBondTypeIndex(species1:Int, species2:Int, structureFileType:Int):Int;
}
