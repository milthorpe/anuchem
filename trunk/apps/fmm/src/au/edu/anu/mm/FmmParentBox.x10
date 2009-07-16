package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a parent node (with one or more child
 * nodes) in the octree division of the 3D simulation space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmParentBox extends FmmBox { 
    /** The non-empty child boxes contained within this box. */
    public val children : GrowableRail[FmmBox];

    public def this(parent : FmmBox) { 
        super(parent);
        this.children = new GrowableRail[FmmBox](0);
    }

    /** Creates a new FmmParentBox with no parent, i.e. a top-level box. */
    public def this() { 
        this(null);
    }
}

