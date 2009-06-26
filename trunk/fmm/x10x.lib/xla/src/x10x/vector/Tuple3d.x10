package x10x.vector;

/** 
 * This class represents a tuple of coordinates.
 * @author milthorpe
 */
public value Tuple3d {
    public val x : double;
    public val y : double;
    public val z : double;

    public def this(x : double, y : double, z : double) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public def toString() : String {
        return ("(" + x + "x + " + y + "y + " + z + "z)");
    }
}

