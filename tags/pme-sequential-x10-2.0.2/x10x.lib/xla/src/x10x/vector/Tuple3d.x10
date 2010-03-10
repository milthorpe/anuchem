package x10x.vector;

/** 
 * This interface represents a 3-tuple of coordinates.
 * @author milthorpe
 */
public abstract class Tuple3d(i : Double, j : Double, k : Double) {
    public def this(i : Double, j : Double, k : Double) {
        property(i, j, k);
    }

    public abstract global safe def toString() : String;
    
    public abstract global def add(b: Tuple3d) : Tuple3d;
    public abstract global def sub(b: Tuple3d) : Tuple3d;
}

