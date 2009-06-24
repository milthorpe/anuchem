package x10x.polar;

import x10x.vector.Point3d;

/** 
 *  This class represents a point in 3D Polar coordinates.
 * @author milthorpe
 */
public value Polar3d {
    public val r: double;
    public val theta: double;
    public val phi: double;

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
        val rxy2 : double = (point.x * point.x) + (point.y * point.y);
        val r2 : double = rxy2 + (point.z * point.z);
        val r : double = Math.sqrt(r2);
        var phi : double;
        var theta : double;
        if (rxy2 == 0.0) {
            if (point.z >= 0.0) {
                theta = 0.0;
            } else {
                theta = Math.PI;
            }
            phi = 0.0;
        } else {
            val rxy : double = Math.sqrt(rxy2);
            phi = Math.atan2(point.y, point.x);
            theta = Math.asin(rxy / r);
        }
        return new Polar3d(r, theta, phi);
    }

    public def toString() : String {
        return ("(r:" + r + ",theta:" + theta + ",phi:" + phi + ")");
    }
}

