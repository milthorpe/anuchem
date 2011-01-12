/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.compiler.CompilerFlags;

/**
 * This class represents a Morton- or Z-ordering of a dense 3-dimensional cubic array.  
 * Each place is assigned a contiguous segment of the Morton curve.
 * NOTE: assumes that the region is cubic (each dimension is the same length)
 * and zero-based.
 */
public final class MortonDist extends Dist(3) {
    /**
     * The place group for this distribution
     */
    private val pg:PlaceGroup;

    /**
     * The number of binary digits per dimension in the Z-index. 
     */
    private val dimDigits : Int;

    /**
     * Cached restricted region for the current place.
     */
    private transient var regionForHere:Region(this.rank);

    public class MortonSubregion extends Region {
        val start : Int;
        val end : Int;
        public def this(start : Int, end : Int) : MortonSubregion{self.rank==3} {
            super(3, false, false);
            this.start = start;
            this.end = end;
        }
        
        public def size() = end-start+1;
        public def totalLength() = MortonDist.this.region.size();
        public def isConvex() = true;
        public def isEmpty() = (end < start);
        public def indexOf(pt:Point):Int {
	        if (pt.rank != 3) return -1;
            return MortonDist.this.getMortonIndex(pt) - start;
        }

        // TODO I hope these aren't used
        public def boundingBox(): Region(rank) {throw new UnsupportedOperationException("boundingBox()");}
        protected  def computeBoundingBox(): Region(rank) {throw new UnsupportedOperationException("computeBoundingBox()");}
        public def min():(Int)=>Int = (i:Int)=> 0;
        public def max():(Int)=>Int = (i:Int)=> MortonDist.this.dimDigits;

        public def contains(that: Region(rank)): boolean {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength() == thatMS.totalLength()
                    && this.start <= thatMS.start && this.end >= thatMS.end) {
                    return true;
                }
            }
            return false;
        }

        public def intersection(that: Region(rank)) : Region(rank) {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength() != thatMS.totalLength()) return Region.makeEmpty(rank); // non-compatible regions
                return MortonDist.this.new MortonSubregion(Math.max(start, thatMS.start), Math.min(end, thatMS.end)) as Region(rank);
            } else {
                throw new UnsupportedOperationException("intersection(Region(!MortonSubregion))");
            }
        }

        public def complement(): Region(rank) {throw new UnsupportedOperationException("complement()");}
        public def product(that:Region): Region {throw new UnsupportedOperationException("product(Region)");}
        public def projection(axis:Int): Region(1) {throw new UnsupportedOperationException("projection(axis : Int)");}
        public def translate(v:Point(rank)): Region(rank){throw new UnsupportedOperationException("translate(Point)");}
        public def eliminate(axis:Int): Region /*(rank-1)*/{throw new UnsupportedOperationException("eliminate(axis : Int)");}

        public def contains(p:Point): boolean {
            if (p.rank != 3) return false;
            val index = MortonDist.this.getMortonIndex(p);
            return (index >= start && index <= end);
        }

        public def iterator(): Iterator[Point(rank)] {
            return new MortonSubregionIterator(this) as Iterator[Point(rank)];
        }

        private class MortonSubregionIterator implements Iterator[Point(3)] {
            var index : Int;

            def this(val r : MortonSubregion) {
                this.index = r.start-1;
            }

            final public def hasNext(): boolean {
                if (index < end) return true;
                else return false;
            }

            final public def next(): Point(3) {
                return MortonDist.this.getPoint(++index);
            } 
        }

        public def scanners():Iterator[Region.Scanner] {
            throw new UnsupportedOperationException("TODO: scanners not defined for MortonSubregion");
        }

        public def toString(): String {
            return "Z[" + start + ".." + end + "]";
        }
    }

    public static def make(r:Region(3)) : MortonDist{self.region==r} {
        return new MortonDist(r, PlaceGroup.WORLD);
    }

    public static def make(r:Region(3), pg:PlaceGroup) : MortonDist{self.region==r} {
        return new MortonDist(r, pg);
    }

    def this(r:Region, pg:PlaceGroup) {
        super(r);
        dimDigits = Math.cbrt(r.size()) as Int;
        this.pg = pg;
    }

    public def places():PlaceGroup = pg;

    public def numPlaces():Int = pg.numPlaces();

    public def get(p:Place):Region(rank) {
        if (p == here) {
            if (regionForHere == null) {
                regionForHere = mortonRegionForPlace(here);
            }
	    return regionForHere;
        } else {
            return mortonRegionForPlace(p);
        }
    }

    private def mortonRegionForPlace(p : Place):Region{self.rank==this.rank} {
        return new MortonSubregion(getPlaceStart(p.id), 
                                   getPlaceEnd(p.id));
    }

    public def regions():Sequence[Region(rank)] {
	    return new Array[Region(rank)](pg.numPlaces(), (i:Int)=>mortonRegionForPlace(pg(i))).sequence();
    }

    public def restriction(r:Region(rank)):Dist(rank) {
	    throw new UnsupportedOperationException("restriction(r:Region(rank))");
    }

    public def restriction(p:Place):Dist(rank) {
	    return Dist.makeConstant(this.get(p));
    }

    public def equals(thatObj:Any):boolean {
	    if (!(thatObj instanceof MortonDist)) return false;
        val that = thatObj as MortonDist;
	    return this.region.size() == that.region.size();
    }

    /**
     * Gets the Morton- or Z-index of a 3-dimensional point.
     * The Morton index is calculated by interleaving binary digits
     * of each dimension.  E.g. 
     * (0, 2, 1) = (00, 10, 01)    = 010001
     * (1, 1, 2) = (01, 01, 10)    = 001110
     * (1, 1, 0) = (01, 01, 00)    = 000110
     */
    public def getMortonIndex(p:Point/*(rank)*/):Int {
        if (p.rank != 3) throw new UnsupportedOperationException("getMortonIndex(p{self.rank!=3})");
        //Console.OUT.println("getMortonIndex for " + p + " dimDigits = " + dimDigits);
        var index : Int = 0;
        var digitMask : Int = Math.pow2(dimDigits-1);
        for (var digit : Int = dimDigits; digit > 0; digit--) {
            for (var dim : Int = 0; dim < 3; dim++) {
                val thisDim = digitMask & p(dim);
                index = index | (thisDim << (digit*2 -dim));
            }
            digitMask = digitMask >> 1;
        }
        //Console.OUT.println(p + " => " + index + " bin = " + index.toBinaryString());
        return index;
    }

    /**
     * Returns a Point(3) corresponding to the given Morton index.
     * @param index the Morton index into the 3D array
     */
    public def getPoint(index:Int) : Point(3) {
        val digitsPerSide = Math.log2(Math.cbrt(region.size()) as Int);
        //Console.OUT.println("getPoint for " + index + " region.size() = " + region.size() + " digitsPerSide = " + digitsPerSide);
        val p = new Array[Int](3);
        var digitMask : Int = region.size() / 2; 
        for (var digit : Int = digitsPerSide; digit > 0; digit--) {
            for (var dim : Int = 0; dim < 3; dim++) {
                val thisDim = digitMask & index;
                p(dim) = p(dim) | (thisDim >> (digit*2 -dim));
                digitMask = digitMask >> 1;
            }
        }
        //Console.OUT.println(index + " => " + p);
        return Point.make(p);
    }

    public def getPlaceStart(placeId : Int) {
        val blockSize = region.size() / Place.MAX_PLACES;
        val numLargerBlocks = region.size() % Place.MAX_PLACES;
        if (placeId < numLargerBlocks) {
            return placeId * (blockSize + 1);
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks) * blockSize;
        }
    }

    public def getPlaceEnd(placeId : Int) {
        val blockSize = region.size() / Place.MAX_PLACES;
        val numLargerBlocks = region.size() % Place.MAX_PLACES;
        if (placeId < numLargerBlocks) {
            return (placeId + 1) * (blockSize + 1) - 1;
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks + 1) * blockSize - 1;
        }
    }

    private final def getPlaceForIndex(index : Int) : Place {
        val blockSize = region.size() / Place.MAX_PLACES;
        val numLargerBlocks = region.size() % Place.MAX_PLACES;
        val firstPart = numLargerBlocks * (blockSize + 1);
        if (index > firstPart) {
            return Place.place(numLargerBlocks + (((index - firstPart) / blockSize) as Int));
        } else {
            return Place.place((index / (blockSize + 1)) as Int);
        }
    } 

    public operator this(pt:Point(rank)):Place {
	    if (CompilerFlags.checkBounds() && !region.contains(pt)) raiseBoundsError(pt);
        val index = getMortonIndex(pt);
        return getPlaceForIndex(index);
    }

    public def offset(pt:Point(rank)):Int {
        // assumes regionForHere is already initialised
        val offset = get(here).indexOf(pt);
        if (offset == -1) {
            if (CompilerFlags.checkBounds() && !region.contains(pt)) raiseBoundsError(pt);
                if (CompilerFlags.checkPlace()) raisePlaceError(pt);
        }
        return offset;
    }

    public def maxOffset() {
        val r = get(here);
        return r.size()-1;
    }

    public def toString(): String {
        var s: String = "MortonDist(";
        var first: boolean = true;
        for (p:Place in places()) {
            if (!first) s += ",";
            s += "" + get(p) + "->" + p.id;
            first = false;
        }
        s += ")";
        return s;
    }

}
