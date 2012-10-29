/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.util.ArrayUtils;
import x10.util.HashMap;

import au.edu.anu.chem.PointCharge;

/**
 * This class represents the Locally Essential Tree (LET) of
 * a single place.  This is the combined interaction lists
 * of all the boxes assigned to that place.
 * @author milthorpe
 */
public class LET {
    public val combinedUList:Rail[UInt];

    /**
     * A cache of multipole copies for the combined V-list of all
     * boxes at this place.  Used to overlap fetching of the multipole
     * expansions with other computation.
     */
    public val multipoleCopies:HashMap[UInt,MultipoleExpansion];

    /**
     * A cache of PointCharge for the combined U-list of all
     * boxes at this place.  Used to store fetched atoms from 
     * non-well-separated boxes for use in direct evaluations 
     * with all atoms at a given place.
     */
    public val cachedAtoms : Rail[Rail[Double]];
    
    public def this(combinedUList:Rail[UInt]) {
        this.combinedUList = combinedUList;
        this.multipoleCopies = new HashMap[UInt,MultipoleExpansion](combinedUList.size); // guess that V-List is at least as big as U-List

        this.cachedAtoms = new Array[Rail[Double]](combinedUList.size);
    }

    public def getMultipoleForOctant(mortonId:UInt) {
        return multipoleCopies.getOrElse(mortonId, null);
    }

    public def setMultipoleForOctant(mortonId:UInt, multipoleExp:MultipoleExpansion) {
        multipoleCopies.put(mortonId, multipoleExp);
    }

    public def getAtomDataForOctant(mortonId:UInt) {
        val cacheIndex = ArrayUtils.binarySearch(combinedUList, mortonId);
        if (cacheIndex >= 0) {
            return cachedAtoms(cacheIndex);
        } else {
            return null;
        }
    }

    public def setAtomDataForOctant(mortonId:UInt, atoms:Rail[Double]) {
        val cacheIndex = ArrayUtils.binarySearch(combinedUList, mortonId);
        if (cacheIndex >= 0) {
            cachedAtoms(cacheIndex) = atoms;
        }
    }
}
