package x10x.vector;

/** 
 * This class represents a tuple of coordinates.
 * @author milthorpe
 */
public value Tuple3d {
    public val i : Double;
    public val j : Double;
    public val k : Double;

    public def this(i : Double, j : Double, k : Double) {
        this.i = i;
        this.j = j;
        this.k = k;
    }

    public def toString() : String {
        return ("(" + i + "i + " + j + "j + " + k + "k)");
    }
    
    public def add(b: Tuple3d) : Tuple3d {
        return new Tuple3d(i + b.i, j + b.j, k + b.k);
    }

    public def sub(b: Tuple3d) : Tuple3d {
        return new Tuple3d(i - b.i, j - b.j, k - b.k);
    }

    public def negate() : Tuple3d {
        return new Tuple3d(-i, -j, -k);
    }
}

