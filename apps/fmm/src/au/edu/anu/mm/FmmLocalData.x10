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

import x10.util.ArrayList;
import x10.util.ArrayUtils;
import x10.util.HashMap;

import au.edu.anu.util.Timer;

public class FmmLocalData {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL:Int = 0;
    public static val TIMER_INDEX_PREFETCH:Int = 1;
    public static val TIMER_INDEX_UPWARD:Int = 2;
    public static val TIMER_INDEX_DOWNWARD:Int = 3;
    public static val TIMER_INDEX_TREE:Int = 4;
    public static val TIMER_INDEX_ASSIGN:Int = 5;
    public static val TIMER_INDEX_REDUCE:Int = 6;
    public static val TIMER_INDEX_REDIST:Int = 7;
    public static val TIMER_INDEX_PARENTS:Int = 8;
    public static val TIMER_INDEX_LET:Int = 9;

    /** The maximum number of levels in the octree. */
    public val dMax:UByte;

    /**
     * The side length of the cube. 
     * To ensure balanced rounding errors within the multipole and local 
     * calculations, all force/energy calculations are performed within 
     * an origin-centred cube.
     */
    public val size:Double; 

    /** All octants held at this place. */
    var octants:HashMap[UInt,Octant];

    /** The load on each leaf octant across all places. */
    val octantLoads:Array[Int];

    /** All leaf octants held at this place. */
    var leafOctants:ArrayList[LeafOctant];

    /** The top-level octants at this place. */
    var topLevelOctants:ArrayList[Octant];

    /** 
     * An array holding the Octant ID of the first leaf octant at each place.
     * The final elements is the ID of the last leaf octant plus one.
     */
    val firstLeafOctant:Rail[UInt] = new Array[UInt](Place.MAX_PLACES+1);

    /** 
     * The locally essential tree at this place. 
     * @see Lashuk et al. (2009).
     */
    var locallyEssentialTree:LET;

    /** The operator arrays for transformations and translations. */
    val fmmOperators:FmmOperators;

    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer:Timer;

    public def this(numTerms:Int, dMax:UByte, ws:Int, size:Double) {
        fmmOperators = new FmmOperators(numTerms, ws);
        this.dMax = dMax;
        this.size = size;
        timer = new Timer(10);
        val maxLeafOctants = Math.pow(8.0, dMax) as Int;
        octantLoads = new Array[Int](maxLeafOctants);
        // TODO construct LET
    }

    /** @return the ID of the place to which the given leaf octant is assigned */
    public def getPlaceId(leafOctantId:OctantId) {
        var placeId:Int = ArrayUtils.binarySearch[UInt](firstLeafOctant, leafOctantId.getLeafMortonId());
        if (placeId < 0) placeId = -(placeId+2);
        return placeId;
    }

    public def getOctant(mortonId:UInt) {
        return octants.getOrElse(mortonId, null);
    }

    public def getOctant(octantId:OctantId) {
        return octants.getOrElse(octantId.getMortonId(), null);
    }

}
