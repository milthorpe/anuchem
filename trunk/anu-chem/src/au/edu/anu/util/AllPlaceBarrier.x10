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
package au.edu.anu.util;

import x10.compiler.Native;

/**
 * This class implements a barrier across all X10 places, 
 * using the X10RT barrier operation.
 * TODO: remove this; this functionality should be provided by clocks,
 * but they are currently inefficient as they are implemented using
 * point-to-point rather than collective communications.
 */
final public class AllPlaceBarrier {
    @Native("c++", "x10rt_barrier()")
    @Native("java", "x10rt_barrier()")
    public static native def barrier() : Int;
}

