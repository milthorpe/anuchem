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
    public val atoms = new GrowableRail[MMAtom]();

    public def this(level : Int, gridLoc : GridLocation, numTerms : Int, parent : FmmBox) { 
        super(level, gridLoc, numTerms, parent);
    }

    public atomic def addAtom(atom : MMAtom) {
        atoms.add(atom);
    }
}

