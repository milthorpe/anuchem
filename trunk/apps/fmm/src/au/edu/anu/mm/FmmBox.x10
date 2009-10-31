package au.edu.anu.mm;

import x10.util.GrowableRail;
import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public global val parent : FmmBox;

    public global val level : Int;

    public global val gridLoc : GridLocation;

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion{self.at(this)};

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion{self.at(this)};

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, gridLoc : GridLocation, numTerms : Int, parent : FmmBox) {
        this.level = level;
        this.gridLoc = gridLoc;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public global def index() : Int {
        dim : Int = Math.pow2(level) as Int;
        return gridLoc.x * dim * dim + gridLoc.y * dim + gridLoc.z;
        //return gridLoc(0) * dim * dim + gridLoc(1) * dim + gridLoc(2);
    }

    public static def getBoxIndex(gridLoc : GridLocation, level : Int) : Int {
        dim : Int = Math.pow2(level);
        return gridLoc.x * dim * dim + gridLoc.y * dim + gridLoc.z;
        //return gridLoc(0) * dim * dim + gridLoc(1) * dim + gridLoc(2);
    }

    public global def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return new Point3d( (gridLoc.x + 0.5) * sideLength - offset,
                            (gridLoc.y + 0.5) * sideLength - offset,
                            (gridLoc.z + 0.5) * sideLength - offset);
        /*
        return new Point3d( (gridLoc(0) + 0.5) * sideLength - offset,
                            (gridLoc(1) + 0.5) * sideLength - offset,
                            (gridLoc(2) + 0.5) * sideLength - offset);
        */
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
        return Math.abs(gridLoc.x - box2GridLoc.x) > ws 
            || Math.abs(gridLoc.y - box2GridLoc.y) > ws 
            || Math.abs(gridLoc.z - box2GridLoc.z) > ws;
        /*
        return Math.abs(gridLoc(0) - box2GridLoc(0)) > ws 
            || Math.abs(gridLoc(1) - box2GridLoc(1)) > ws 
            || Math.abs(gridLoc(2) - box2GridLoc(2)) > ws;
        */
    }

    public global def getTranslationIndex(box2 : FmmBox) : GridLocation {
        return GridLocation(gridLoc.x - box2.gridLoc.x, gridLoc.y - box2.gridLoc.y, gridLoc.z - box2.gridLoc.z);
        //return [gridLoc(0) - box2.gridLoc(0), gridLoc(1) - box2.gridLoc(1), gridLoc(2) - box2.gridLoc(2)];
    }
}

