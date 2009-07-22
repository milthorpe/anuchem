package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmParentBox extends FmmBox {    
    public val children : GrowableRail[FmmBox];

    public def this(level : Int, location : ValRail[Int]{length==3},numTerms : Int, parent : FmmBox) { 
        super(level, location, numTerms, parent);
        children = new GrowableRail[FmmBox](0);
    }
}

