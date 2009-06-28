package x10x.vector;

/** 
 * This class represents a 3D vector
 * @author V.Ganesh
 */
public value Vector3d extends Tuple3d {
    public static val NULL = new Vector3d(0, 0, 0);

    public static val UX = new Vector3d(1.0, 0, 0);
    public static val UY = new Vector3d(0, 1.0, 0);
    public static val UZ = new Vector3d(0, 0, 1.0);
    
    public def this(x : Double, y : Double, z : Double) {
        super(x, y, z);
    }

    public def dot(vec : Vector3d) : Double {
        return this.x * vec.x + this.y * vec.y + this.z * vec.z;
    }

    public def cross(vec : Vector3d) : Vector3d {
        return new Vector3d(y * vec.z - z * vec.y,
                            z * vec.x - x * vec.z,
                            x * vec.y - y * vec.x);
    }

    public def mul(k : Double) : Vector3d {
        return new Vector3d(this.x * k, this.y * k, this.z * k);
    }

    public def mixedProduct(v2 : Vector3d, v3 : Vector3d) : Double {
        return this.dot(v2.cross(v3));
    }

    public def lengthSquared() : Double {
        return x * x + y * y + z * z;
    }

    public def length() : Double {
        return Math.sqrt(x * x + y * y + z * z);
    }

    /* Should which name to use - length or magnitude? */
    public def magnitude() : Double {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public def angleWith(vec : Vector3d) : Double {
        val aDotb = this.dot(vec);
        val ab    = magnitude() * vec.magnitude();
        
        return Math.acos(aDotb / ab);
    }

    public def normalize() : Vector3d {
        val norm = 1.0 / length();

        return new Vector3d(x * norm, y * norm, z * norm);
    }

    public def negate() : Vector3d {
        return new Vector3d(-x, -y, -z);
    }
}

