package x10x.vector;

/** 
 * This interface represents a 3-tuple of coordinates.
 * @author milthorpe
 */
public interface Tuple3d(i : Double, j : Double, k : Double) {
    public global safe def toString() : String;
    
    public global def add(b: Tuple3d) : Tuple3d;
    public global def sub(b: Tuple3d) : Tuple3d;
}

