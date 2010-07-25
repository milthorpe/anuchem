/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

public class ExpansionRegion extends Region  {
    // XTENLANG-49
    static type ExpansionRegion(rank:Int) = ExpansionRegion{self.rank==rank};
    global val p : Int;

    /**
     * Constructs a new ExpansionRegion for a multipole
     * or local expansion of the given dimension.
     * @param p the dimension of the expansion
     */
    public def this(p : Int): ExpansionRegion(2) {
        super(2, false, true);
        this.p = p;
    }

    public global def isConvex() : Boolean {
        return true;
    }

    public global def isEmpty() : Boolean {
        return false;
    }

    public global safe def min(): ValRail[int] = [0,-p];

    public global safe def max(): ValRail[int] = [p,p];

    /**
     * Returns the number of points in this region.
     * These expansions have a peculiarly shaped region(abs(x1)<=x0 && 0<=x0<=p)
     * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3),
     * which gives a size of (p+1)^2.
     */
    public global def size() : Int { 
        return (p+1) * (p+1);
    }

    public global def contains(p: Point): boolean {
        if (p.rank == 2) {
            return (p(0) >= 0 && p(0) <= this.p && Math.abs(p(1)) <= p(0));
        }
        throw U.unsupported(this, "contains(" + p + ")");
    }

    public global def contains(r:Region(rank)): boolean {
        // TODO
        throw U.unsupported(this, "contains(Region)");
    }

    public global def complement(): Region(rank) {
        // TODO
        throw U.unsupported(this, "complement()");
    }

    public global def intersection(t: Region(rank)): Region(rank) {
        // TODO
        throw U.unsupported(this, "intersection()");
    }

    public global def product(r: Region): Region {
        // TODO
        throw U.unsupported(this, "product()");
    }

    public global def translate(v: Point(rank)): Region(rank) {
        // TODO
        throw U.unsupported(this, "translate()");
    }

    public global def projection(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return 0..p;
            case 1:
                return -p..p;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public global def eliminate(axis: int): Region(1) {
        switch (axis) {
            case 0:
                return -p..p;
            case 1:
                return 0..p;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public global def boundingBox(): Region(rank) {
        return [0..p,-p..p];
    }

    protected global def computeBoundingBox(): Region(rank) {
        return [0..p,-p..p];
    }

    public global def iterator(): Iterator[Point(rank)] {
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

        final public safe def hasNext(): boolean {
            if (l <= p && m <= l) return true;
            else return false;
        }

        final public safe def next(): Point(2) {
            nextPoint : Point(2) = Point.make(l,m);
            if (m < l) m++;
            else {
                l++;
                m = -l; 
            }
            return nextPoint;
        } 
    }

    public global def scanners():Iterator[Region.Scanner]! {
        throw U.unsupported("TODO: scanners not defined for ExpansionRegion");
    }

    public global safe def toString(): String {
        return "ExpansionRegion (p = " + p + ")";
    }
}
