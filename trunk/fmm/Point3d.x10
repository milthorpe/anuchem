/**
 * @author Josh Milthorpe
 */

package au.edu.anu.math;

/** 
 *  This class represents a point in 3D cartesian space.
 */
public value Point3d {
    val x : double;
    val y: double;
    val z : double;

    public def this(x : double, y : double, z : double) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

