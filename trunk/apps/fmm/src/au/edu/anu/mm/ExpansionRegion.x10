package au.edu.anu.mm;

import x10.array.BaseRegion;
import x10.array.U;

public value class ExpansionRegion extends BaseRegion{rank==2} {
    val p : Int;

    /**
     * Constructs a new ExpansionRegion for a multipole
     * or local expansion of the given dimension.
     * @param p the dimension of the expansion
     */
    public def this(p : Int){p>0} {
        super(p, false, true);
        this.p = p;
    }

    public def isConvex() : Boolean {
        return true;
    }

    public def isEmpty() : Boolean {
        return false;
    }

    public def min(): ValRail[int] {
        return [0,-p];
    }

    public def max(): ValRail[int] {
        return [p,p];
    }   

    /**
     * Returns the number of points in this region.
     * These expansions have a peculiarly shaped region(abs(x0)<=x1 && x0<=p)
     * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3),
     * which gives a size of p^2 + 2(p+1).
     */
    public def size() : Int { 
        return p * p + 2 * (p + 1);
    }

    public def contains(p: Point): boolean {
        if (p.rank == 2) {
            return (p(0) >= 0 && p(0) <= this.p && Math.abs(p(1)) <= p(0));
        }
        throw U.unsupported(this, "contains(" + p + ")");
    }

    public def complement(): Region(rank) {
        // TODO
        throw U.unsupported(this, "complement()");
    }

    public def intersection(t: Region(rank)): Region(rank) {
        // TODO
        throw U.unsupported(this, "intersection()");
    }

    public def product(r: Region): Region {
        // TODO
        throw U.unsupported(this, "product()");
    }

    public def projection(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return 0..p;
            case 1:
                return -p..p;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public def eliminate(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return -p..p;
            case 1:
                return 0..p;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public def boundingBox(): Region(rank) {
        return [0..p,-p..p];
    }

    public def computeBoundingBox(): Region(rank) {
        return [0..p,-p..p];
    }

    public def iterator(): Iterator[Point(rank)] {
        return new ExpansionRegionIterator(this) as Iterator[Point(rank)];
    }

    private class ExpansionRegionIterator implements Iterator[Point(2)] {
        val p : Int;
        var l : Int;
        var m : Int;

        def this(val r : ExpansionRegion) {
            this.p = r.p;
            this.l = 0;
            this.m = 0;
        }

        final public def hasNext(): boolean {
            if (l <= p && m <= l) return true;
            else return false;
        }

        final public def next(): Point(2) {
            nextPoint : Point(2) = [l,m];
            if (m < l) m++;
            else {
                l++;
                m = -l; 
            }
            return nextPoint;
        } 
    }

}
