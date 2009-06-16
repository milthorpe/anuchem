/**
 * @author Josh Milthorpe
 */

package x10x.vector;

/** 
 *  This class represents a point in 3D cartesian space.
 */
public value Point3d {
    public val x : double;
    public val y: double;
    public val z : double;

    public def this(x : double, y : double, z : double) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

