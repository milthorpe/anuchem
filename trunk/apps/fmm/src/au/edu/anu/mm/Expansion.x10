package au.edu.anu.mm;

import x10.lang.Complex;
import x10.util.StringBuilder;

/**
 * This is the superclass for multipole and local expansions, as used in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * These expansions have a peculiarly shaped region(abs(x0)<=x1 && x0<=p)
 * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3), which is constructed
 * here by subtracting two halfspaces from a rectangular region. 
 * @author milthorpe
 */
public value Expansion {
    /** The terms X_{lm} (with m >= 0) in this expansion */
    public val terms : Array[Complex]{rank==2};

    public def this(p : Int) {
        var expRegion : Region{rank==2} = [0..p,-p..p];
        expRegion = expRegion - Region.makeHalfspace([1,1],1);
        expRegion = expRegion - Region.makeHalfspace([1,-1],1);
        this.terms = Array.make[Complex](expRegion->here, (Point)=> Complex.ZERO);
    }

    public def add(e : Expansion) {
        val p : Int = terms.region.max(0);
        for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-i..i]) {
                this.terms(i,j) = this.terms(i,j).add(e.terms(i,j));
            }
        }
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
