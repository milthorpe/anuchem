/**
 * @author Josh Milthorpe
 */

package x10x.polar;

import x10x.vector.Point3d;

/** 
 *  This class represents a point in 3D Polar coordinates.
 */
public value Polar3d {
    val r: double;
    val theta: double;
    val phi: double;

    public def this(r : double, theta : double, phi : double) {
        this.r = r;
        this.theta = theta;
        this.phi = phi;
    }

    /* Returns a cartesian representation of this point. */
    public def toPoint3d() {
        sineTheta : double = Math.sin(theta);
        return new Point3d(r * sineTheta * Math.cos(phi), 
                            r * sineTheta * Math.sin(phi),
                            r * Math.cos(theta)
                           );
    }

    /** Returns a polar representation of the given cartesian point. */
    public static def getPolar3d(point : Point3d) {
        r : double = Math.sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
        return new Polar3d(r,
                           Math.atan(point.y / point.x),
                           Math.acos(point.z / r)
                          );
    } 
}

