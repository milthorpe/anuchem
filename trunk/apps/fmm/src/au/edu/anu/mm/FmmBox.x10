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

    public static global def getBoxLocation(index : Int, level : Int) : GridLocation {
        dim : Int = Math.pow2(level) as Int;
        return GridLocation(index / (dim * dim), (index / dim) % dim, index % dim);
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
        return wellSeparated(ws, box2.gridLoc);
    }

    /**
     * Returns true if this box is well-separated from <code>box2Loc</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global def wellSeparated(ws : Int, box2Loc : GridLocation) : Boolean {
        if (level < 2)
            return false;
        // TODO can do reduction on a Rail?
        return Math.abs(gridLoc.x - box2Loc.x) > ws 
            || Math.abs(gridLoc.y - box2Loc.y) > ws 
            || Math.abs(gridLoc.z - box2Loc.z) > ws;
    }

    public global def getTranslationIndex(box2 : FmmBox) : GridLocation {
        return getTranslationIndex(box2.gridLoc);
    }

    public global def getTranslationIndex(loc2 : GridLocation) : GridLocation {
        return GridLocation(gridLoc.x - loc2.x, gridLoc.y - loc2.y, gridLoc.z - loc2.z);
    }
}

