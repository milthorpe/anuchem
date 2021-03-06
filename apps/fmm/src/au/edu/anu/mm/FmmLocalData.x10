/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012-2013.
 */

package au.edu.anu.mm;

import x10.util.ArrayList;
import x10.util.RailUtils;
import x10.util.HashMap;
import x10.util.HashSet;

import au.edu.anu.util.Timer;
import au.edu.anu.chem.mm.ParticleData;

public class FmmLocalData {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL = 0;
    public static val TIMER_INDEX_PREFETCH = 1;
    public static val TIMER_INDEX_UPWARD = 2;
    public static val TIMER_INDEX_M2L = 3;
    public static val TIMER_INDEX_DOWNWARD = 4;
    public static val TIMER_INDEX_TREE = 5;
    public static val TIMER_INDEX_SORT = 6;
    public static val TIMER_INDEX_BALANCE = 7;
    public static val TIMER_INDEX_REDIST = 8;
    public static val TIMER_INDEX_PARENTS = 9;
    public static val TIMER_INDEX_LET = 10;

    /** The maximum number of levels in the octree. */
    public var dMax:UByte;

    /**
     * The side length of the cube. 
     * To ensure balanced rounding errors within the multipole and local 
     * calculations, all force/energy calculations are performed within 
     * an origin-centred cube.
     */
    public var size:Double;

    /**
     * Particle data: positions, charges and mass of particles at this place.
     */
    public var particleData:ParticleData;

    /** All octants held at this place. */
    val octants = new HashMap[UInt,Octant]();

    /** The load on each leaf octant across all places. */
    var octantLoads:Rail[Long];

    /** All leaf octants held at this place. */
    val leafOctants = new ArrayList[LeafOctant]();

    /** The top-level octants at this place. */
    val topLevelOctants = new ArrayList[Octant]();

    /** 
     * A Rail holding the Octant ID of the first leaf octant at each place.
     * The final elements is the ID of the last leaf octant plus one.
     */
    val firstLeafOctant:Rail[UInt] = new Rail[UInt](Place.numPlaces()+1);

    /** 
     * The locally essential tree at this place. 
     * @see Lashuk et al. (2009).
     */
    var locallyEssentialTree:LET;

    /** The operator arrays for transformations and translations. */
    var fmmOperators:FmmOperators;

    /** cost estimates (in ns) per interaction = [P2P, M2L] */
    val cost:Rail[Long] = new Rail[Long](2);
    public static val ESTIMATE_P2P=0;
    public static val ESTIMATE_M2L=1;

    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(11);

    public def init(numTerms:Int, dMax:UByte, ws:Int, size:Double) {
        this.dMax = dMax;
        this.size = size;
        fmmOperators = new FmmOperators(numTerms, ws);
        val maxLeafOctants = Math.pow(8.0, dMax as Int) as Int;
        octantLoads = new Rail[Long](maxLeafOctants);
        // TODO construct LET
    }

    /** @return the ID of the place to which the given leaf octant is assigned */
    public def getPlaceId(leafOctantId:OctantId) {
        return getPlaceId(leafOctantId.getMortonId());
    }

    /** @return the ID of the place to which the given leaf octant is assigned */
    public def getPlaceId(mortonId:UInt) {
        var placeId:Long = RailUtils.binarySearch[UInt](firstLeafOctant, mortonId & OctantId.LEAF_MASK);
        if (placeId < 0) placeId = -(placeId+2);
        return placeId;
    }

    public def getOctant(mortonId:UInt) {
        return octants.getOrElse(mortonId, null);
    }

    public def getOctant(octantId:OctantId) {
        return octants.getOrElse(octantId.getMortonId(), null);
    }

    public def getCombinedUList(ws:Int) {
        val combinedUSet = new HashSet[UInt]();
        for(octant in leafOctants) {
            val uList = octant.getUList();
            for (mortonId in uList) {
                combinedUSet.add(mortonId);
            }
        }

        //Console.OUT.println("at " + here + " combined U-list:");
        val combinedUList = new Rail[UInt](combinedUSet.size());
        var j:Long = 0;
        for (mortonId in combinedUSet) {
            combinedUList(j++) = mortonId;
            //Console.OUT.println(mortonId);
        }
        RailUtils.sort(combinedUList);
        return combinedUList;
    }
}
