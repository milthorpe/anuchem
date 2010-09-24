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

import x10.util.GrowableRail;

import x10x.vector.Point3d;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    public val atoms : GrowableRail[MMAtom] = new GrowableRail[MMAtom]();

    /** The U-list consists of all leaf boxes not well-separated to this box. */
    private var uList : ValRail[Point(3)];

    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : GlobalRef[FmmBox]) { 
        super(level, x, y, z, numTerms, parent);
    }

    public atomic def addAtom(atom : MMAtom) {
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

