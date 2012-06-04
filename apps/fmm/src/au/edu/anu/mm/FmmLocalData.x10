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

import au.edu.anu.util.Timer;

public class FmmLocalData {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL:Int = 0;
    public static val TIMER_INDEX_PREFETCH:Int = 1;
    public static val TIMER_INDEX_UPWARD:Int = 2;
    public static val TIMER_INDEX_DOWNWARD:Int = 3;
    public static val TIMER_INDEX_TREE:Int = 4;
    public static val TIMER_INDEX_PLACEHOLDER:Int = 5;

    /** All leaf octants held at this place. */
    var leafOctants:ArrayList[LeafOctant];

    /** The top-level shared octant at this place. */
    var parentOctant:SharedOctant;

    /** 
     * The locally essential tree at this place. 
     * @see Lashuk et al. (2009).
     */
    var locallyEssentialTree:LET;

    /** The operator arrays for transformations and translations. */
    val fmmOperators:FmmOperators;

    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer:Timer;

    public def this(numTerms:Int, ws:Int) {
        fmmOperators = new FmmOperators(numTerms, ws);
        timer = new Timer(6);
        // TODO construct LET
    }

}
