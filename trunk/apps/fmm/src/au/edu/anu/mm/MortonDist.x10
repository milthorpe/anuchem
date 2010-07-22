package au.edu.anu.mm;

import x10.array.BaseDist;
import x10.array.BaseRegion;
import x10.array.U;

/**
 * This class represents a Morton- or Z-ordering of a dense 3-dimensional cubic array.  
 * Each place is assigned a contiguous segment of the Morton curve.
 * NOTE: assumes that the region is cubic (each dimension is the same length)
 * and zero-based.
 */
public class MortonDist extends BaseDist{self.rank==3} {
    global val totalLength : Int;

    public static class MortonSubregion extends BaseRegion {
        global val start : Int;
        global val end : Int;
        global val totalLength : Int;
        public def this(start : Int, end : Int, totalLength : Int) : MortonSubregion{self.rank==3} {
            super(3, false, false);
            this.start = start;
            this.end = end;
            this.totalLength = totalLength;
        }
        
        public global def size() = end-start+1;
        public global def isConvex() = true;
        public global def isEmpty() = (end < start);
        // TODO I hope these aren't used
        public global def boundingBox(): Region(rank) {throw U.unsupported(this, "boundingBox()");}
        global protected  def computeBoundingBox(): Region(rank) {throw U.unsupported(this, "computeBoundingBox()");}
        public global def min(): ValRail[int] = ValRail.make(rank, (Int) => 0);
        public global def max(): ValRail[int] = ValRail.make(rank, (Int) => (Math.cbrt(totalLength) as Int));

        public global def contains(that: Region(rank)): boolean {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength == thatMS.totalLength
                    && this.start <= thatMS.start && this.end >= thatMS.end) {
                    return true;
                }
            }
            return false;
        }

        public global def intersection(that: Region(rank)) : Region(rank) {
            if (that instanceof MortonSubregion) {
                val thatMS = that as MortonSubregion;
                if (this.totalLength != thatMS.totalLength) return Region.makeEmpty(rank); // non-compatible regions
                return new MortonSubregion(Math.max(start, thatMS.start), Math.min(end, thatMS.end), this.totalLength) as Region(rank);
            } else {
                throw U.unsupported(this, "intersection(Region(!MortonSubregion))");
            }
        }

        public global def complement(): Region(rank) {throw U.unsupported(this, "complement()");}
        public global def product(that: Region): Region { throw U.unsupported(this, "product(Region)");}
        public global def projection(axis: int): Region(1) {throw U.unsupported(this, "projection(axis : Int)");}
        public global def translate(v: Point(rank)): Region(rank){throw U.unsupported(this, "translate(Point)");}
        public global def eliminate(axis: int): Region /*(rank-1)*/{throw U.unsupported(this, "eliminate(axis : Int)");}

        public global def contains(p: Point): boolean {
            if (p.rank != 3) return false;
            val index = MortonDist.getMortonIndex(p, this.totalLength);
            return (index >= start && index <= end);
        }

        public global def iterator(): Iterator[Point(rank)] {
            return new MortonSubregionIterator(this) as Iterator[Point(rank)];
        }

        private class MortonSubregionIterator implements Iterator[Point(3)] {
            val totalLength : Int;
            var index : Int;

            def this(val r : MortonSubregion) {
                this.totalLength = r.totalLength;
                this.index = r.start-1;
            }

            final public safe def hasNext(): boolean {
                if (index < end) return true;
                else return false;
            }

            final public safe def next(): Point(3) {
                return MortonDist.getPoint(++index, this.totalLength);
            } 
        }

        public global safe def toString(): String {
            return "Z[" + start + ".." + end + "]";
        }
    }

    public static def make(r: Region(3)) : MortonDist{self.region==r} {
        return MortonDist.make(r, Place.places);
    }

    public static def make(r: Region(3), ps: ValRail[Place]) : MortonDist{self.region==r} {
        val totalLength = r.size();
        val init = (p:Int) => new MortonSubregion(getPlaceStart(p,ps.length(),totalLength), 
                                                  getPlaceEnd(p,ps.length(),totalLength), 
                                                  totalLength) as Region(3);
        val subregions = ValRail.make[Region(3)](ps.length(), init);
        return new MortonDist(r, ps, subregions);
    }

    public def this(r: Region, ps: ValRail[Place], rs: ValRail[Region(r.rank)]): MortonDist{self.region==r} {
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
    public static safe def getMortonIndex(p : Point/*(rank)*/, totalLength : Int) : Int {
        if (p.rank != 3) throw U.unsupported("getMortonIndex(p{self.rank!=3})");
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
    public static safe def getPoint(index : Int, totalLength : Int) : Point(3) {
        val digitsPerSide = Math.log2(Math.cbrt(totalLength) as Int);
        //Console.OUT.println("getPoint for " + index + " totalLength = " + totalLength + " digitsPerSide = " + digitsPerSide);
        val p = Rail.make[Int](3);
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

    public static safe def getPlaceStart(placeId : Int, numPlaces : Int, totalLength : Int) {
        val blockSize = totalLength / numPlaces;
        val numLargerBlocks = totalLength % numPlaces;
        if (placeId < numLargerBlocks) {
            return placeId * (blockSize + 1);
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks) * blockSize;
        }
    }

    public static safe def getPlaceEnd(placeId : Int, numPlaces : Int, totalLength : Int) {
        val blockSize = totalLength / numPlaces;
        val numLargerBlocks = totalLength % numPlaces;
        if (placeId < numLargerBlocks) {
            return (placeId + 1) * (blockSize + 1) - 1;
        } else {
            val firstPortion = numLargerBlocks * (blockSize + 1);
            return firstPortion + (placeId - numLargerBlocks + 1) * blockSize - 1;
        }
    }

    public global safe def apply(pt: Point/*(rank)*/): Place {
        if (pt.rank != 3) throw U.unsupported("getMortonIndex(p{self.rank!=3})");
        val index = getMortonIndex(pt, totalLength);
        for (p:Place in places) {
            val mr = get(p) as MortonSubregion;
            if (mr.contains(pt)) {
                return p;
            }
        }
        throw new ArrayIndexOutOfBoundsException("point " + pt + " not contained in distribution");
    }

    public global safe def toString(): String {
        var s: String = "MortonDist(";
        var first: boolean = true;
        for (p:Place in places) {
            if (!first) s += ",";
            s +=  get(p) + "->" + p.id;
            first = false;
        }
        s += ")";
        return s;
    }

}
