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

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents the Locally Essential Tree (LET) of
 * a single place.  This is the combined interaction lists
 * of all the boxes assigned to that place.
 * @author milthorpe
 */
public class LocallyEssentialTree {
    public val combinedUList : ValRail[Point(3)];
    public val combinedVList : ValRail[ValRail[Point(3)]];
    public val uListMin : ValRail[Int](3);
    public val uListMax : ValRail[Int](3);
    public val vListMin : ValRail[ValRail[Int](3)];
    public val vListMax : ValRail[ValRail[Int](3)];

    /**
     * A cache of multipole copies for the combined V-list of all
     * boxes at this place.  Used to overlap fetching of the multipole
     * expansions with other computation.
     * The ValRail has one element for each level; each element
     * holds the portion of the combined V-list for that level.
     */
    public global val multipoleCopies : ValRail[PeriodicArray[MultipoleExpansion](3)!];

    /**
     * A cache of packed for the combined U-list of all
     * boxes at this place.  Used to store fetched packed atoms
     * from non-well-separated boxes for use in direct evaluations 
     * with all atoms at a given place.
     * @see FmmLeafBox.getPackedAtoms()
     */
    public global val packedAtoms : PeriodicArray[ValRail[MMAtom.PackedRepresentation]](3)!;
    
    public def this(combinedUList : ValRail[Point(3)],
                combinedVList : ValRail[ValRail[Point(3)]],
                uListMin : ValRail[Int](3),
                uListMax : ValRail[Int](3),
                vListMin : ValRail[ValRail[Int](3)],
                vListMax : ValRail[ValRail[Int](3)]) {
        this.combinedUList = combinedUList;
        this.combinedVList = combinedVList;
        this.uListMin = uListMin;
        this.uListMax = uListMax;
        this.vListMin = vListMin;
        this.vListMax = vListMax;
        val multipoleCopies = Rail.make[PeriodicArray[MultipoleExpansion](3)!](combinedVList.length());
        for ([i] in 0..combinedVList.length()-1) {
            if (combinedVList(i) != null) {
                val multipoleCopiesLevelRegion : Region(3) = [vListMin(i)(0)..vListMax(i)(0), vListMin(i)(1)..vListMax(i)(1), vListMin(i)(2)..vListMax(i)(2)];
                multipoleCopies(i) = new PeriodicArray[MultipoleExpansion](multipoleCopiesLevelRegion);
            }
        }
        this.multipoleCopies = ValRail.make(multipoleCopies);

        val packedAtomsRegion : Region(3) = [uListMin(0)..uListMax(0), uListMin(1)..uListMax(1), uListMin(2)..uListMax(2)];
        this.packedAtoms = new PeriodicArray[ValRail[MMAtom.PackedRepresentation]](packedAtomsRegion);
    }
}
