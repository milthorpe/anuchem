package x10x.vector;

/** 
 * This class represents a 3D vector
 * @author V.Ganesh
 */
public class Vector3d extends Tuple3d {
    public static val NULL = new Vector3d(0, 0, 0);

    public static val UX = new Vector3d(1.0, 0, 0);
    public static val UY = new Vector3d(0, 1.0, 0);
    public static val UZ = new Vector3d(0, 0, 1.0);
    
    public def this(i : Double, j : Double, k : Double) {
        super(i, j, k);
    }

    public global def dot(vec : Vector3d) : Double {
        return this.i * vec.i + this.j * vec.j + this.k * vec.k;
    }

    public global def cross(vec : Vector3d) : Vector3d {
        return new Vector3d(j * vec.k - k * vec.j,
                            k * vec.i - i * vec.k,
                            i * vec.j - j * vec.i);
    }

    public global def mul(c : Double) : Vector3d {
        return new Vector3d(this.i * c, this.j * c, this.k * c);
    }

    public global def mixedProduct(v2 : Vector3d, v3 : Vector3d) : Double {
        return this.dot(v2.cross(v3));
    }

    public global def lengthSquared() : Double {
        return i * i + j * j + k * k;
    }

    public global def length() : Double {
        return Math.sqrt(i * i + j * j + k * k);
    }

    /* Should which name to use - length or magnitude? */
    public global def magnitude() : Double {
        return Math.sqrt(i * i + j * j + k * k);
    }

    public global def angleWith(vec : Vector3d) : Double {
        val aDotb = this.dot(vec);
        val ab    = magnitude() * vec.magnitude();
        
        return Math.acos(aDotb / ab);
    }

    public global def normalize() : Vector3d {
        val norm = 1.0 / length();

        return new Vector3d(i * norm, j * norm, k * norm);
    }

    public global def negate() : Vector3d {
        return new Vector3d(-i, -j, -k);
    }
}

