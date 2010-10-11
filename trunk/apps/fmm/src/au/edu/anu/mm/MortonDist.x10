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

import x10.array.BaseDist;

/**
 * This class represents a Morton- or Z-ordering of a dense 3-dimensional cubic array.  
 * Each place is assigned a contiguous segment of the Morton curve.
 * NOTE: assumes that the region is cubic (each dimension is the same length)
 * and zero-based.
 */
public class MortonDist extends BaseDist{self.rank==3} {
    val totalLength : Int;

    public static class MortonSubregion extends Region {
        val start : Int;
        val end : Int;
        val totalLength : Int;
        public def this(start : Int, end : Int, totalLength : Int) : MortonSubregion{self.rank==3} {
            super(3, false, false);
            this.start = start;
            this.end = end;
            this.totalLength = totalLength;
        }
        
        public def size() = end-start+1;
        public def isConvex() = true;
        public def isEmpty() = (end < start);
        public def indexOf(pt:Point) {
	        if (pt.rank != 3) return -1;
            return MortonDist.getMortonIndex(pt, this.totalLength);
        }

        // TODO I hope these aren't used
        public def boundingBox(): Region(rank) {throw new UnsupportedOperationException("boundingBox()");}
        protected  def computeBoundingBox(): Region(rank) {throw new UnsupportedOperationException("computeBoundingBox()");}
        public def min():(int)=>int = (i:int)=> 0;
        public def max():(int)=>int = (i:int)=> (Math.cbrt(totalLength) as Int);

        public def contains(that: Region(rank)): boolean {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength == thatMS.totalLength
                    && this.start <= thatMS.start && this.end >= thatMS.end) {
                    return true;
                }
            }
            return false;
        }

        public def intersection(that: Region(rank)) : Region(rank) {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength != thatMS.totalLength) return Region.makeEmpty(rank); // non-compatible regions
                return new MortonSubregion(Math.max(start, thatMS.start), Math.min(end, thatMS.end), this.totalLength) as Region(rank);
            } else {
                throw new UnsupportedOperationException("intersection(Region(!MortonSubregion))");
            }
        }

        public def complement(): Region(rank) {throw new UnsupportedOperationException("complement()");}
        public def product(that: Region): Region {throw new UnsupportedOperationException("product(Region)");}
        public def projection(axis: int): Region(1) {throw new UnsupportedOperationException("projection(axis : Int)");}
        public def translate(v: Point(rank)): Region(rank){throw new UnsupportedOperationException("translate(Point)");}
        public def eliminate(axis: int): Region /*(rank-1)*/{throw new UnsupportedOperationException("eliminate(axis : Int)");}

        public def contains(p: Point): boolean {
            if (p.rank != 3) return false;
            val index = MortonDist.getMortonIndex(p, this.totalLength);
            return (index >= start && index <= end);
        }

        public def iterator(): Iterator[Point(rank)] {
            return new MortonSubregionIterator(this) as Iterator[Point(rank)];
        }

        private class MortonSubregionIterator implements Iterator[Point(3)] {
            val totalLength : Int;
            var index : Int;

            def this(val r : MortonSubregion) {
                this.totalLength = r.totalLength;
                this.index = r.start-1;
            }

            final public def hasNext(): boolean {
                if (index < end) return true;
                else return false;
            }

            final public def next(): Point(3) {
                return MortonDist.getPoint(++index, this.totalLength);
            } 
        }

        public def scanners():Iterator[Region.Scanner] {
            throw new UnsupportedOperationException("TODO: scanners not defined for MortonSubregion");
        }

        public def toString(): String {
            return "Z[" + start + ".." + end + "]";
        }
    }

    public static def make(r: Region(3)) : MortonDist{self.region==r} {
        return MortonDist.make(r, new Array[Place](Place.MAX_PLACES, (i : Int) => Place.place(i)));
    }

    public static def make(r: Region(3), ps: Array[Place]{rail}) : MortonDist{self.region==r} {
        val totalLength = r.size();
        val init = (p:Int) => new MortonSubregion(getPlaceStart(p,ps.size,totalLength), 
                                                  getPlaceEnd(p,ps.size,totalLength), 
                                                  totalLength) as Region(3);
        val subregions = new Array[Region(3)](ps.size, init);
        return new MortonDist(r, ps, subregions);
    }

    public def this(r: Region, ps: Array[Place](1), rs: Array[Region(r.rank)](1)): MortonDist{self.region==r} {
        super(r, ps, rs);
        totalLength = r.size();
    }

    /**
     * Gets the Morton- or Z-index of a 3-dimensional point.
     * The Morton index is calculated by interleaving binary digits
     * of each dimension.  E.g. 
     * (0, 2, 1) = (00, 10, 01)    = 010001
     * (1, 1, 2) = (01, 01, 10)    = 001110
     * (1, 1, 0) = (01, 01, 00)    = 000110
     */
    public static def getMortonIndex(p : Point/*(rank)*/, totalLength : Int) : Int {
        if (p.rank != 3) throw new UnsupportedOperationException("getMortonIndex(p{self.rank!=3})");
        val digitsPerSide = Math.cbrt(totalLength) as Int;
        //Console.OUT.println("getMortonIndex for " + p + " digitsPerSide = " + digitsPerSide);
        var index : Int = 0;
        var digitMask : Int = Math.pow2(digitsPerSide-1);
        for (var digit : Int = digitsPerSide; digit > 0; digit--) {
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
     * @param totalLength the total length of the array = (2n)^3
     */
    public static def getPoint(index : Int, totalLength : Int) : Point(3) {
        val digitsPerSide = Math.log2(Math.cbrt(totalLength) as Int);
        //Console.OUT.println("getPoint for " + index + " totalLength = " + totalLength + " digitsPerSide = " + digitsPerSide);
        val p = new Array[Int](3);
        var digitMask : Int = totalLength / 2; 
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

    public static def getPlaceStart(placeId : Int, numPlaces : Int, totalLength : Int) {
        val blockSize = totalLength / numPlaces;
        val numLargerBlocks = totalLength % numPlaces;
        if (placeId < numLargerBlocks) {
            return placeId * (blockSize + 1);
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks) * blockSize;
        }
    }

    public static def getPlaceEnd(placeId : Int, numPlaces : Int, totalLength : Int) {
        val blockSize = totalLength / numPlaces;
        val numLargerBlocks = totalLength % numPlaces;
        if (placeId < numLargerBlocks) {
            return (placeId + 1) * (blockSize + 1) - 1;
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks + 1) * blockSize - 1;
        }
    }

    public def apply(pt: Point/*(rank)*/): Place {
        if (pt.rank != 3) throw new UnsupportedOperationException("getMortonIndex(p{self.rank!=3})");
        val index = getMortonIndex(pt, totalLength);
        for (p:Place in places()) {
            val mr = get(p) as MortonSubregion;
            if (mr.contains(pt)) {
                return p;
            }
        }
        throw new ArrayIndexOutOfBoundsException("point " + pt + " not contained in distribution");
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
