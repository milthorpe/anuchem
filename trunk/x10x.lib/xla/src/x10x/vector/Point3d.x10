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
        this(t.i(), t.j(), t.k());
    }

    public global safe def toString() = ("(" + i + "i + " + j + "j + " + k + "k)");

    /**
     * @return a new point displaced from this point by the vector <code>b</code>
     */
    public global safe def add(b: Vector3d) : Point3d {
        return Point3d(i + b.i(), j + b.j(), k + b.k());
    }

    public safe operator this + (that:Vector3d):Point3d {
        return this.add(that);
    }

    /**
     * @return the vector from that point to this point
     */
    public global safe def vector(b: Point3d) : Vector3d {
        return Vector3d(i - b.i(), j - b.j(), k - b.k());
    }   

    public safe operator this - (that:Point3d):Vector3d {
        return this.vector(that);
    }

    /**
     * check equality of two points
     */
    public static safe operator (x:Point3d) == (y:Point3d) = ((x.i==y.i) && (x.j==y.j) && (x.k==y.k));

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

