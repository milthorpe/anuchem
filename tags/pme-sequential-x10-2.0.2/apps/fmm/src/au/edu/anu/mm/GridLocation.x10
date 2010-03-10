package au.edu.anu.mm;

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

