package x10x.vector;

/** 
 * This class represents a tuple of coordinates.
 * @author milthorpe
 */
public class Tuple3d {
    public global val i : Double;
    public global val j : Double;
    public global val k : Double;

    public def this(i : Double, j : Double, k : Double) {
        this.i = i;
        this.j = j;
        this.k = k;
    }

    public global safe def toString() : String {
        return ("(" + i + "i + " + j + "j + " + k + "k)");
    }
    
    public global def add(b: Tuple3d) : Tuple3d {
        return new Tuple3d(i + b.i, j + b.j, k + b.k);
    }

    public global def sub(b: Tuple3d) : Tuple3d {
        return new Tuple3d(i - b.i, j - b.j, k - b.k);
    }

    public global def negate() : Tuple3d {
        return new Tuple3d(-i, -j, -k);
    }
}

