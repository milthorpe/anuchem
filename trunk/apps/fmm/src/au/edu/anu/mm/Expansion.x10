package au.edu.anu.mm;

import x10.lang.Complex;
import x10x.vector.Tuple3d;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.polar.Polar3d;

/**
 * This is the superclass for multipole and local expansions, as used in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public value Expansion {
    /** The terms X_{lm} (with m >= 0) in this expansion */
    public val terms : Array[Complex]{rank==2};

    public def this(p : Int) {
        this.terms = Array.make[Complex]([0..p, -p..p]->here, (val (i,j): Point)=> Complex.ZERO);
    }
}
