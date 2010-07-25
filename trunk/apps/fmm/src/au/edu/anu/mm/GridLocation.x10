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
package au.edu.anu.mm;

/**
 * This class represents a box location in the three dimensional
 * grid.  Boxes are numbered 0..dim-1 in each dimension.
 */
public struct GridLocation {
    public val x : Int;
    public val y : Int;
    public val z : Int;

    public def this(x : Int, y : Int, z : Int) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

