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

import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public global val parent : FmmBox;

    public global val level : Int;
    public global val x : Int;
    public global val y : Int;
    public global val z : Int;

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private var vList : ValRail[Point(3)];

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion{self.at(this)};

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion{self.at(this)};

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : FmmBox) {
        this.level = level;
        this.x = x;
        this.y = y;
        this.z = z;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public global safe def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return Point3d( (x + 0.5) * sideLength - offset,
                        (y + 0.5) * sideLength - offset,
                        (z + 0.5) * sideLength - offset);
    }

    /**
     * Returns true if this box is well-separated from <code>x,y,z</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global safe def wellSeparated(ws : Int, x2 : Int, y2 : Int, z2 : Int) : Boolean {
        return Math.abs(x - x2) > ws 
            || Math.abs(y - y2) > ws 
            || Math.abs(z - z2) > ws;
    }

    /**
     * Returns true if this box is well-separated from <code>box2</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global safe def wellSeparated(ws : Int, box2 : FmmBox) : Boolean {
        return Math.abs(x - box2.x) > ws 
            || Math.abs(y - box2.y) > ws 
            || Math.abs(z - box2.z) > ws;
    }

    public global safe def getTranslationIndex(level2 : Int, x2 : Int, y2 : Int, z2 : Int) : Point(4) {
        return Point.make(level, x2-x, y2-y, z2-z);
    }

    public global safe def getTranslationIndex(box2 : FmmBox) : Point(4) {
        return Point.make(level, box2.x-x, box2.y-y, box2.z-z);
    }

    public def getVList() = this.vList;

    public def setVList(vList : ValRail[Point(3)]) {
        this.vList = vList;
    }

    public global safe def toString(): String {
        return "FmmBox level " + level + " (" + x + "," + y + "," + z + ")";
    }
}

