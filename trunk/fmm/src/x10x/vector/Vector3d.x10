package x10x.vector;

/** 
 * This class represents a 3D vector
 * (Initial DRAFT)
 *
 * @author V.Ganesh
 */
public value Vector3d {
    public val i : Double;
    public val j : Double;
    public val k : Double;

    public static val NULL = new Vector3d(0, 0, 0);

    public static val UX = new Vector3d(1.0, 0, 0);
    public static val UY = new Vector3d(0, 1.0, 0);
    public static val UZ = new Vector3d(0, 0, 1.0);
    
    public def this(i : Double, j : Double, k : Double) {
        this.i = i;
        this.j = j;
        this.k = k;
    }

    public def dot(vec:Vector3d) : Double {
        var result:Double = 0.0;
                
        result += this.i * vec.i;
        result += this.j * vec.j;
        result += this.k * vec.k;
        
        return result;
    }

    public def cross(vec:Vector3d) : Vector3d {
        return new Vector3d(j * vec.k - k * vec.j,
                            k * vec.i - i * vec.k,
                            i * vec.j - j * vec.i);
    }

    public def mul(k:Double) : Vector3d {
        return new Vector3d(this.i * k, this.j * k, this.k * k);
    }

    public def mixedProduct(v2:Vector3d, v3:Vector3d) : Double {
        return this.dot(v2.cross(v3));
    }

    public def magnitude() : Double {
        var length:Double = 0.0;
        
        length += i * i;
        length += j * j;
        length += k * k;
        
        return Math.sqrt(length);
    }

    public def angleWith(vec:Vector3d) : Double {
        val aDotb = this.dot(vec);
        val ab    = magnitude() * vec.magnitude();
        
        return Math.acos(aDotb / ab);
    }

    public def normalize() : Vector3d {
        val mag = magnitude();

        return new Vector3d(i / mag, j / mag, k / mag);
    }

    public def negate() : Vector3d {
        return new Vector3d(-i, -j, -k);
    }

    public def toString() : String {
        return ("(" + i + "i + " + j + "j + " + k + "k)");
    }
}

