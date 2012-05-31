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

import x10.util.HashMap;

import au.edu.anu.chem.PointCharge;

/**
 * This class represents the Locally Essential Tree (LET) of
 * a single place.  This is the combined interaction lists
 * of all the boxes assigned to that place.
 * @author milthorpe
 */
public class LET {
    public val combinedUList:Rail[OctantId];
    public val combinedVList:Rail[OctantId];

    /**
     * A cache of multipole copies for the combined V-list of all
     * boxes at this place.  Used to overlap fetching of the multipole
     * expansions with other computation.
     */
    public val multipoleCopies:HashMap[OctantId,MultipoleExpansion];

    /**
     * A cache of PointCharge for the combined U-list of all
     * boxes at this place.  Used to store fetched atoms from 
     * non-well-separated boxes for use in direct evaluations 
     * with all atoms at a given place.
     */
    public val cachedAtoms : Rail[Rail[PointCharge]];
    
    public def this(combinedUList:Rail[OctantId],
                combinedVList:Rail[OctantId]) {
        this.combinedUList = combinedUList;
        this.combinedVList = combinedVList;
        this.multipoleCopies = new HashMap[OctantId,MultipoleExpansion](combinedVList.size);

        this.cachedAtoms = new Array[Rail[PointCharge]](combinedUList.size);
    }
}
