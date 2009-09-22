package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    public val atoms = new GrowableRail[Atom]();

    public def this(level : Int, location : ValRail[Int]{length==3},numTerms : Int, parent : FmmBox) { 
        super(level, location, numTerms, parent);
    }

    public def addAtom(atom : Atom) {
        atoms.add(atom);
    }
}

