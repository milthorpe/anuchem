/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
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

