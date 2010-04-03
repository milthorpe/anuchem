package au.edu.anu.mm;

import x10x.vector.Point3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public global val parent : FmmBox;

    public global val level : Int;
    public global val x : Int;
    public global val y : Int;
    public global val z : Int;

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private var vList : ValRail[Point(3)];

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion{self.at(this)};

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion{self.at(this)};

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : FmmBox) {
        this.level = level;
        this.x = x;
        this.y = y;
        this.z = z;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public global safe def index() : Int {
        dim : Int = Math.pow2(level) as Int;
        return x * dim * dim + y * dim + z;
    }

    public global safe def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return new Point3d( (x + 0.5) * sideLength - offset,
                            (y + 0.5) * sideLength - offset,
                            (z + 0.5) * sideLength - offset);
    }

    /**
     * Returns true if this box is well-separated from <code>x,y,z</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global safe def wellSeparated(ws : Int, x2 : Int, y2 : Int, z2 : Int) : Boolean {
        if (level < 2)
            return false;
        return Math.abs(x - x2) > ws 
            || Math.abs(y - y2) > ws 
            || Math.abs(z - z2) > ws;
    }

    /**
     * Returns true if this box is well-separated from <code>boxIndex</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global safe def wellSeparated(ws : Int, boxIndex : Point(3)) : Boolean {
        if (level < 2)
            return false;
        return Math.abs(x - boxIndex(0)) > ws 
            || Math.abs(y - boxIndex(1)) > ws 
            || Math.abs(z - boxIndex(2)) > ws;
    }

    /**
     * Returns true if this box is well-separated from <code>box2</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public global safe def wellSeparated(ws : Int, box2 : FmmBox) : Boolean {
        if (level < 2)
            return false;
        return Math.abs(x - box2.x) > ws 
            || Math.abs(y - box2.y) > ws 
            || Math.abs(z - box2.z) > ws;
    }

    public global safe def getTranslationIndex(boxIndex2 : Point(4)) : Point(4) {
        return Point.make(level, x - boxIndex2(1), y - boxIndex2(2), z - boxIndex2(3));
    }

    public global safe def getTranslationIndex(level2 : Int, x2 : Int, y2 : Int, z2 : Int) : Point(4) {
        return Point.make(level, x - x2, y - y2, z - z2);
    }

    public global safe def getTranslationIndex(box2 : FmmBox) : Point(4) {
        return Point.make(level, x - box2.x, y - box2.y, z - box2.z);
    }

    /**
     * TODO this should not be necessary once XTENLANG-787 is resolved
     * @return a local copy at the current place of this box's multipole expansion
     */
    public global def getMultipoleExpansionLocalCopy(p : Int) : MultipoleExpansion! {
        val data = at (this) {Expansion.getData(p, multipoleExp)};
        return new MultipoleExpansion(p, data);
    }

    /**
     * TODO this should not be necessary once XTENLANG-787 is resolved
     * @return a local copy at the current place of this box's local expansion
     */
    public global def getLocalExpansionLocalCopy(p : Int) : LocalExpansion! {
        val data = at (this) {Expansion.getData(p, localExp)};
        return new LocalExpansion(p, data);
    }

    public def getVList() = this.vList;

    public def setVList(vList : ValRail[Point(3)]) {
        this.vList = vList;
    }

    public global safe def toString(): String {
        return "FmmBox level " + level + " (" + x + "," + y + "," + z + ")";
    }
}

