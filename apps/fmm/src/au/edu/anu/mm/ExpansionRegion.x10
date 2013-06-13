/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.mm;

import x10.regionarray.Region;

public class ExpansionRegion extends Region{rect} {
// TODO not really rect! should be 'dense' XTENLANG-3000
    static type ExpansionRegion(rank:Int) = ExpansionRegion{self.rank==rank,rect==true};
    val p:Long;

    /**
     * Constructs a new ExpansionRegion for a multipole
     * or local expansion of the given dimension.
     * @param p the dimension of the expansion
     */
    public def this(p : Int): ExpansionRegion(2) {
        super(2, true, false);
        this.p = p;
    }

    public def isConvex() : Boolean {
        return true;
    }

    public def isEmpty() : Boolean {
        return false;
    }

    public def min():(Int)=>Long = (i:Int)=> {
        if (i==0) return 0L;
        else if (i==1) return -p;
        else throw new ArrayIndexOutOfBoundsException("min: "+i+" is not a valid rank for "+this);
    };

    public def max():(Int)=>Long = (i:Int)=> {
        if (i==0) return p;
        else if (i==1) return p;
        else throw new ArrayIndexOutOfBoundsException("max: "+i+" is not a valid rank for "+this);
    };

    /**
     * Returns the number of points in this region.
     * These expansions have a peculiarly shaped region(abs(x1)<=x0 && 0<=x0<=p)
     * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3),
     * which gives a size of (p+1)^2.
     */
    public def size():Long { 
        return (p+1) * (p+1);
    }

    public def contains(p: Point): boolean {
        if (p.rank == 2) {
            return (p(0) >= 0 && p(0) <= this.p && Math.abs(p(1)) <= p(0));
        }
        throw new UnsupportedOperationException("contains(" + p + ")");
    }

    public def contains(r:Region(rank)): boolean {
        if (r instanceof ExpansionRegion && (r as ExpansionRegion).p == this.p)
            return true;
        throw new UnsupportedOperationException("contains(Region)");
    }

    public def complement(): Region(rank) {
        // TODO
        throw new UnsupportedOperationException("complement()");
    }

    public def intersection(t: Region(rank)): Region(rank) {
        // TODO
        throw new UnsupportedOperationException("intersection()");
    }

    public def product(that: Region): Region{self!=null} {
        // TODO
        throw new UnsupportedOperationException("product()");
    }

    public def translate(v: Point(rank)): Region(rank) {
        // TODO
        throw new UnsupportedOperationException("translate()");
    }

    public def projection(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return Region.make(0, p);
            case 1:
                return Region.make(-p, p);
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def eliminate(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return Region.make(-p, p);
            case 1:
                return Region.make(0, p);
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def indexOf(pt:Point) {
	    if (pt.rank != 2) return -1L;
        return indexOf(pt(0), pt(1));
    }

    public final def indexOf(i0:Long, i1:Long) {
        return i0*(i0+1) + i1;
    }

    public def boundingBox(): Region(rank) {
        return computeBoundingBox();
    }

    protected def computeBoundingBox(): Region(rank) {
        val r = Region.make(0..p, -p..p);
        return r;
    }

    public def iterator(): Iterator[Point(rank)] {
        return new ExpansionRegionIterator(this) as Iterator[Point(rank)];
    }

    private class ExpansionRegionIterator implements Iterator[Point(2)] {
        val p:Long;
        var l:Long;
        var m:Long;

        def this(val r : ExpansionRegion) {
            this.p = r.p;
            this.l = 0L;
            this.m = 0L;
        }

        final public def hasNext(): boolean {
            if (l <= p && m <= l) return true;
            else return false;
        }

        final public def next(): Point(2) {
            nextPoint : Point(2) = Point.make(l,m);
            if (m < l) m++;
            else {
                l++;
                m = -l; 
            }
            return nextPoint;
        } 
    }

    public def toString(): String {
        return "ExpansionRegion (p = " + p + ")";
    }
}
