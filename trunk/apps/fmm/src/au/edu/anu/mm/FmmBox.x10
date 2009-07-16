package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public abstract class FmmBox { 
    /** The multipole expansion of the charges within this box. */
    public var multipoleExp : MultipoleExpansion;

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public var localExp : LocalExpansion;

    /** The parent box containing this box. */
    public val parent : FmmBox;

    public def this(parent : FmmBox) { 
        this.parent = parent;
    }
}

