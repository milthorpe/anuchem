package x10x.vector;

import x10.compiler.Inline;

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

    property i() = i;
    property j() = j;
    property k() = k;

    public def toString() = ("(" + i + "i + " + j + "j + " + k + "k)");

    /**
     * @return a new point displaced from this point by the vector <code>b</code>
     */
    public @Inline def add(b: Vector3d) : Point3d {
        return Point3d(i + b.i(), j + b.j(), k + b.k());
    }

    public @Inline operator this + (that:Vector3d):Point3d {
        return this.add(that);
    }

    /**
     * @return this point scaled by the given factor
     */
    public @Inline operator this * (that:Double):Point3d {
        return this.scale(that);
    }

    public @Inline def scale(c : Double) : Point3d {
        return Point3d(this.i * c, this.j * c, this.k * c);
    }

    /**
     * Performs anisotropic (non-uniform) scaling of this point.
     * @return this point scaled by the vector v
     */
    public @Inline def scale(v : Vector3d) : Point3d {
        return Point3d(this.i * v.i, this.j * v.j, this.k * v.k);
    }

    /**
     * @return the vector from that point to this point
     */
    public @Inline def vector(b: Point3d) : Vector3d {
        return Vector3d(i - b.i, j - b.j, k - b.k);
    }   

    public @Inline operator this - (that:Point3d):Vector3d {
        return this.vector(that);
    }

    public @Inline def distanceSquared(b : Point3d) {
        val di = i - b.i;
        val dj = j - b.j;
        val dk = k - b.k;
        return di*di
             + dj*dj
             + dk*dk;
    }  

    public @Inline def distance(b : Point3d) {
        return Math.sqrt(distanceSquared(b));
    }
}

