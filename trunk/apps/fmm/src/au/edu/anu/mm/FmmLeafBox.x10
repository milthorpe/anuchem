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

import x10.util.GrowableRail;

import x10x.vector.Point3d;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    public val atoms : GrowableRail[MMAtom{self.at(this)}]{self.at(this)} = new GrowableRail[MMAtom{self.at(this)}]();

    /** The U-list consists of all leaf boxes not well-separated to this box. */
    private var uList : ValRail[Point(3)];

    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : FmmBox) { 
        super(level, x, y, z, numTerms, parent);
    }

    public safe atomic def addAtom(atom : MMAtom{self.at(this)}) {
        atoms.add(atom);
    }

    public def getUList() = this.uList;

    public def setUList(uList : ValRail[Point(3)]) {
        this.uList = uList;
    }
    
    /*
     * Returns atom charges and coordinates in packed representation
     */
    public def getPackedAtoms() : ValRail[MMAtom.PackedRepresentation] {
        if (atoms.length() > 0) {
            return ValRail.make[MMAtom.PackedRepresentation](atoms.length(), (i : Int) => atoms(i).getPackedRepresentation());
        } else {
            return null;
        }
    }
}

