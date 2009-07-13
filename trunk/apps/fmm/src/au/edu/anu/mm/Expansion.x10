package au.edu.anu.mm;

import x10.lang.Complex;
import x10.util.StringBuilder;

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

    public def toString() : String {
        val p : Int = terms.region.max(0);
        s : StringBuilder = new StringBuilder();
        for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-i..i]) {
			    s.add(terms(i,j) + " ");
            }
            s.add("\n");
		}
        return s.toString();
    }
}
