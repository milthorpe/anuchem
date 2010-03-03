package x10x.vector;

/** 
 * This class represents a point in 3D cartesian space.
 * @author milthorpe
 */
public struct Point3d(i : Double, j : Double, k : Double) implements Tuple3d {
    public def this(i : Double, j : Double, k : Double) {
        property(i, j, k);
    }

    public def this(t : Tuple3d) {
        property(t.i(), t.j(), t.k());
    }

    public global safe def toString() = ("(" + i + "i + " + j + "j + " + k + "k)");
    
    public global safe def add(b: Tuple3d) : Tuple3d {
        return Point3d(i + b.i(), j + b.j(), k + b.k());
    }

    public global safe def sub(b: Tuple3d) : Tuple3d {
        return Point3d(i - b.i(), j - b.j(), k - b.k());
    }

    // TODO - specific sub-typed methods should not be required
    public global safe def sub(b: Point3d) : Point3d {
        return Point3d(i - b.i(), j - b.j(), k - b.k());
    }

    public global safe def distanceSquared(b : Point3d) {
        return (i - b.i) * (i - b.i)
             + (j - b.j) * (j - b.j)
             + (k - b.k) * (k - b.k);
    }  

    public global safe def distance(b : Point3d) {
        return Math.sqrt((i - b.i) * (i - b.i)
                       + (j - b.j) * (j - b.j)
                       + (k - b.k) * (k - b.k));
    }
}

