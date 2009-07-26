package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public val parent : FmmBox;

    public val level : Int;

    public val location : ValRail[Int]{length==3};

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion;

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion;

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, location : ValRail[Int]{length==3}, numTerms : Int, parent : FmmBox) {
        this.level = level;
        this.location = location;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

        public def index() : Int {
        dim : Int = Math.pow2(level) as Int;
        return location(0) * dim * dim + location(1) * dim + location(2);
    }

    public static def getBoxIndex(location : ValRail[Int]{length==3}, level : Int) : Int {
        dim : Int = Math.pow2(level);
        return location(0) * dim * dim + location(1) * dim + location(2);
    }

    public def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return new Point3d( (location(0) + 0.5) * sideLength - offset,
                            (location(1) + 0.5) * sideLength - offset,
                            (location(2) + 0.5) * sideLength - offset);
    }

    /**
     * Returns true if this box is well-separated from <code>box2</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    def wellSeparated(ws : Int, box2 : FmmBox) : Boolean {
        if (level < 2)
            return false;
        if (this == box2)
            return false;
        // TODO can do reduction on a Rail?
        return Math.abs(location(0) - box2.location(0)) > ws 
            || Math.abs(location(1) - box2.location(1)) > ws 
            || Math.abs(location(2) - box2.location(2)) > ws;
    }
}

