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

import x10.util.StringBuilder;

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
    public val terms : Array[Complex](2)!;

    public def this(p : Int) {
        //var expRegion : Region(2) = [0..p,-p..p];
        //expRegion = expRegion - Region.makeHalfspace([1,1],1);
        //expRegion = expRegion - Region.makeHalfspace([1,-1],1);
        val expRegion = new ExpansionRegion(p);
        this.terms = new Array[Complex](expRegion, (Point) => Complex.ZERO);
    }

    public def this(p : Int, data : ValRail[Complex]) {
        val expRegion = new ExpansionRegion(p);
	// FIXME: DG.  Change this constructor to take an Array to copyFrom to avoid useless zeroing in new LocalRectArray.
	    val tmp = new Array[Complex](expRegion);
        //data.copyTo(0, tmp.raw(), 0, data.length());
        for ((i,j) in tmp) {
            tmp(i,j) = data(i*i + i+j);
        }
        this.terms = tmp;
    }

    public atomic def add(e : Expansion!) {
        for ((i,j) in terms.region) {
            this.terms(i,j) = this.terms(i,j) + e.terms(i,j);
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
        val size = (p+1)*(p+1);
        val data = Rail.make[Complex](size);
        for ((i,j) in source.terms) {
            data(i*i + i+j) = source.terms(i,j);
        }
        return data;
        //return source.terms.raw();
    }
}
