package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public global val parent : Box[FmmBox];

    public global val level : Int;

    public global val gridLoc : ValRail[Int]{length==3};

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion{self.at(this)};

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion{self.at(this)};

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, gridLoc : ValRail[Int]{length==3}, numTerms : Int, parent : Box[FmmBox]) {
        this.level = level;
        this.gridLoc = gridLoc;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public global def index() : Int {
        dim : Int = Math.pow2(level) as Int;
        return gridLoc(0) * dim * dim + gridLoc(1) * dim + gridLoc(2);
    }

    public static def getBoxIndex(gridLoc : ValRail[Int]{length==3}, level : Int) : Int {
        dim : Int = Math.pow2(level);
        return gridLoc(0) * dim * dim + gridLoc(1) * dim + gridLoc(2);
    }

    public global def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return new Point3d( (gridLoc(0) + 0.5) * sideLength - offset,
                            (gridLoc(1) + 0.5) * sideLength - offset,
                            (gridLoc(2) + 0.5) * sideLength - offset);
    }

    /**
     * Returns true if this box is well-separated from <code>box2</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global def wellSeparated(ws : Int, box2 : FmmBox) : Boolean {
        if (level < 2)
            return false;
        //if (this == box2)
        //    return false;
        // TODO can do reduction on a Rail?
        val box2GridLoc = box2.gridLoc;
        return Math.abs(gridLoc(0) - box2GridLoc(0)) > ws 
            || Math.abs(gridLoc(1) - box2GridLoc(1)) > ws 
            || Math.abs(gridLoc(2) - box2GridLoc(2)) > ws;
    }

    public global def getTranslationIndex(box2 : FmmBox) : ValRail[Int]{length==3} {
        val box2GridLoc = box2.gridLoc;
        return [gridLoc(0) - box2GridLoc(0), gridLoc(1) - box2GridLoc(1), gridLoc(2) - box2GridLoc(2)];
    }
}

