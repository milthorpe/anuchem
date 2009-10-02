package x10x.vector;

/** 
 * This class represents a point in 3D cartesian space.
 * @author milthorpe
 */
public value Point3d extends Tuple3d {
    public def this(i : Double, j : Double, k : Double) {
        super(i, j, k);
    }

    public def distanceSquared(b : Point3d) {
        return (i - b.i) * (i - b.i)
             + (j - b.j) * (j - b.j)
             + (k - b.k) * (k - b.k);
    }  

    public def distance(b : Point3d) {
        return Math.sqrt((i - b.i) * (i - b.i)
                       + (j - b.j) * (j - b.j)
                       + (k - b.k) * (k - b.k));
    }
}

