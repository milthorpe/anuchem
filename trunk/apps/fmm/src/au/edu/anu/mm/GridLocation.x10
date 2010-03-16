package au.edu.anu.mm;

/**
 * This class represents a box location in the three dimensional
 * grid.  Boxes are numbered 0..dim-1 in each dimension.
 */
public struct GridLocation {
    public val x : Int;
    public val y : Int;
    public val z : Int;

    public def this(x : Int, y : Int, z : Int) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

