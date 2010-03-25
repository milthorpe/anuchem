package au.edu.anu.mm;

import x10.util.StringBuilder;
import x10.array.FastArray;

/**
 * This is the superclass for multipole and local expansions, as used in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * These expansions have a peculiarly shaped region(abs(x0)<=x1 && 0<=x1<=p)
 * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3), which is constructed
 * here by subtracting two halfspaces from a rectangular region. 
 * @author milthorpe
 */
public class Expansion {
    /** The terms X_{lm} (with m >= 0) in this expansion */
    public val terms : Array[Complex](2);

    public def this(p : Int) {
        //var expRegion : Region(2) = [0..p,-p..p];
        //expRegion = expRegion - Region.makeHalfspace([1,1],1);
        //expRegion = expRegion - Region.makeHalfspace([1,-1],1);
        val expRegion = new ExpansionRegion(p);
        this.terms = Array.make[Complex](expRegion->here, (Point) => Complex.ZERO);
    }

    public def this(p : Int, data : ValRail[Complex]) {
        val expRegion = new ExpansionRegion(p);
        this.terms = Array.make[Complex](expRegion->here);
        data.copyTo(0, (this.terms as FastArray[Complex]).raw(), 0, data.length());
    }

    public atomic def add(e : Expansion!) {
        for (val p in terms.region) {
            this.terms(p) = this.terms(p) + e.terms(p);
        }
    }

    public global safe def toString() : String {
        return at(this) {toStringLocal()};
    }

    private safe def toStringLocal() : String {
        val p : Int = terms.region.max(0);
        val s = new StringBuilder();
        for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-i..i]) {
			    s.add(terms(i,j) + " ");
            }
            s.add("\n");
		}
        return s.toString();
    }

    /**
     * TODO this should not be necessary once XTENLANG-787 is resolved
     * @return the expansion terms, shoehorned into a ValRail
     */
    public static safe def getData(p : Int, source : Expansion!) : ValRail[Complex] {
        // HACK - use raw data - this won't work without NO_CHECKS!!
        return (source.terms as FastArray[Complex]).raw();
    }

    /** The nearest integer square root below <code>n</code> */
    protected static safe def squareRootFloor(n : Int) = (Math.sqrt(n as Double) as Int);
}
