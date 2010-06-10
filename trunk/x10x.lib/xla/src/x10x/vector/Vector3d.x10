package x10x.vector;

/** 
 * This class represents a vector in 3D cartesian space.
 * @author V.Ganesh, milthorpe
 */
public struct Vector3d(i : Double, j : Double, k : Double) implements Tuple3d {
    public static val NULL = Vector3d(0, 0, 0);

    public static val UX = Vector3d(1.0, 0, 0);
    public static val UY = Vector3d(0, 1.0, 0);
    public static val UZ = Vector3d(0, 0, 1.0);
    
    public def this(i : Double, j : Double, k : Double) {
        property(i, j, k);
    }

    public def this(t : Tuple3d) {
        this(t.i(), t.j(), t.k());
    }

    public global safe def toString() = ("(" + i + "i + " + j + "j + " + k + "k)");

    /**
     * @return the sum of this vector and the given vector
     */
    public safe operator this + (that:Vector3d):Vector3d {
        return this.add(that);
    }

    public global safe def add(b: Vector3d) : Vector3d {
        return Vector3d(i + b.i(), j + b.j(), k + b.k());
    }

    /**
     * @return the difference of vector x and vector y
     */
    public static safe operator (x:Vector3d) - (y:Vector3d):Vector3d {
        return x.sub(y);
    }

    public global safe def sub(b: Vector3d) : Vector3d {
        return Vector3d(i - b.i(), j - b.j(), k - b.k());
    }

    /**
     * @return the dot product of this vector and the given vector
     */
    public safe operator this * (that:Vector3d):Double {
        return this.dot(that);
    }

    public global safe def dot(vec : Vector3d) : Double {
        return this.i * vec.i + this.j * vec.j + this.k * vec.k;
    }

    public global safe def cross(vec : Vector3d) : Vector3d {
        return Vector3d(j * vec.k - k * vec.j,
                        k * vec.i - i * vec.k,
                        i * vec.j - j * vec.i);
    }

    /**
     * @return the product of this vector and the given scalar
     */
    public safe operator this * (that:Double):Vector3d {
        return this.mul(that);
    }

    /**
     * @return the product of the given double and the given vector
     */
    public static safe operator (x:Double) * (y:Vector3d): Vector3d = y * x;

    /**
     * check equality of two vectors    
     */
    public static safe operator (x:Vector3d) == (y:Vector3d) = ((x.i==y.i) && (x.j==y.j) && (x.k==y.k));

    public global safe def mul(c : Double) : Vector3d {
        return Vector3d(this.i * c, this.j * c, this.k * c);
    }

    public global safe def mixedProduct(v2 : Vector3d, v3 : Vector3d) : Double {
        return this.dot(v2.cross(v3));
    }

    public global safe def lengthSquared() : Double {
        return i * i + j * j + k * k;
    }

    public global safe def length() : Double {
        return Math.sqrt(i * i + j * j + k * k);
    }

    public global safe def maxNorm() : Double {
        return Math.max(Math.max(Math.abs(i), Math.abs(j)), Math.abs(k));
    }

    /* Should which name to use - length or magnitude? */
    public global safe def magnitude() : Double {
        return Math.sqrt(i * i + j * j + k * k);
    }

    public global safe def angleWith(vec : Vector3d) : Double {
        val aDotb = this.dot(vec);
        val ab    = magnitude() * vec.magnitude();
        
        return Math.acos(aDotb / ab);
    }

    public global def normalize() : Vector3d {
        val norm = 1.0 / length();
        return Vector3d(i * norm, j * norm, k * norm);
    }

    /** 
     * Returns the geometric inverse of this vector a^(-1) = a / ||a||^2
     */
    public global safe def inverse() : Vector3d {
        val l2 = lengthSquared();
        return Vector3d(i / l2, j / l2, k / l2);
    }

    public static safe operator - (x:Vector3d) = x.negate();

    public global safe def negate() : Vector3d {
        return Vector3d(-i, -j, -k);
    }
}

