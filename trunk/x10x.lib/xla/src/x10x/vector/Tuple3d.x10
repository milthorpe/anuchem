package x10x.vector;

/** 
 * This class represents a tuple of coordinates.
 * @author milthorpe
 */
public value Tuple3d {
    public val x : Double;
    public val y : Double;
    public val z : Double;

    public def this(x : Double, y : Double, z : Double) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public def toString() : String {
        return ("(" + x + "x + " + y + "y + " + z + "z)");
    }
}

