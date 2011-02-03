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
    public val combinedUList : Array[Point(3)](1){rect,rail};
    public val combinedVList : Array[Array[Point(3)](1){rect,rail}](1){rect,rail};
    public val uListMin : Array[Int](1){rect,rail};
    public val uListMax : Array[Int](1){rect,rail};
    public val vListMin : Array[Array[Int](1){rect,rail}](1){rect,rail};
    public val vListMax : Array[Array[Int](1){rect,rail}](1){rect,rail};

    /**
     * A cache of multipole copies for the combined V-list of all
     * boxes at this place.  Used to overlap fetching of the multipole
     * expansions with other computation.
     * The Array has one element for each level; each element
     * holds the portion of the combined V-list for that level.
     */
    public val multipoleCopies : Array[DistArray[MultipoleExpansion](3)](1){rect,rail};

    /**
     * A cache of packed for the combined U-list of all
     * boxes at this place.  Used to store fetched packed atoms
     * from non-well-separated boxes for use in direct evaluations 
     * with all atoms at a given place.
     * @see FmmLeafBox.getPackedAtoms()
     */
    public val packedAtoms : DistArray[Array[MMAtom.PackedRepresentation](1){rect,rail}](3);
    
    public def this(combinedUList : Array[Point(3)](1){rect,rail},
                combinedVList : Array[Array[Point(3)](1){rect,rail}](1){rect,rail},
                uListMin : Array[Int](1){rect,rail},
                uListMax : Array[Int](1){rect,rail},
                vListMin : Array[Array[Int](1){rect,rail}](1){rect,rail},
                vListMax : Array[Array[Int](1){rect,rail}](1){rect,rail}) {
        this.combinedUList = combinedUList;
        this.combinedVList = combinedVList;
        this.uListMin = uListMin;
        this.uListMax = uListMax;
        this.vListMin = vListMin;
        this.vListMax = vListMax;
        val multipoleCopies = new Array[DistArray[MultipoleExpansion](3)](combinedVList.size);
        for ([i] in 0..(combinedVList.size-1)) {
            if (combinedVList(i) != null) {
                val multipoleCopiesLevelRegion = vListMin(i)(0)..vListMax(i)(0) * vListMin(i)(1)..vListMax(i)(1) * vListMin(i)(2)..vListMax(i)(2);
                multipoleCopies(i) = DistArray.make[MultipoleExpansion](new PeriodicDist(Dist.makeConstant(multipoleCopiesLevelRegion)));
            }
        }
        this.multipoleCopies = multipoleCopies;

        val packedAtomsRegion = uListMin(0)..uListMax(0) * uListMin(1)..uListMax(1) * uListMin(2)..uListMax(2);
        this.packedAtoms = DistArray.make[Array[MMAtom.PackedRepresentation](1){rect,rail}](new PeriodicDist(Dist.makeConstant(packedAtomsRegion)));
    }
}
