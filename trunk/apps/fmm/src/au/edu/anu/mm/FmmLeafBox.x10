package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    public val atoms : GrowableRail[Atom]{self.home==this.home} = new GrowableRail[Atom]();

    public def this(level : Int, gridLoc : GridLocation, numTerms : Int, parent : FmmBox) { 
        super(level, gridLoc, numTerms, parent);
    }

    public atomic def addAtom(atom : Atom) {
        atoms.add(atom);
    }
}

