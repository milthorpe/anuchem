package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a leaf node in the 3D octree division of 
 * space for the fast multipole method.  There are no child boxes.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox { 
    /**
     * Creates a new FmmLeafBox.
     */
    public def this(parent : FmmBox) { 
        super(parent);
    }
}

