package x10.regionarray;

public class TriangularRegion extends Region{rect} {
// TODO not really rect! should be 'dense' XTENLANG-3000
    // XTENLANG-49
    static type TriangularRegion(rank:Int) = TriangularRegion{self.rank==rank,rect};

    val dim:Long;
    val rowMin:Long;
    val colMin:Long;
    val lower:Boolean;

    /**
     * Constructs a new TriangularRegion 
     * @param p the dimension of the triangular array
     */
    public def this(rowMin:Long, colMin:Long, size:Long, lower: Boolean): TriangularRegion(2) {
        super(2, true, (rowMin==0L&&colMin==0L));
        this.dim = size;
        this.rowMin = rowMin;
        this.colMin = colMin;
        this.lower = lower;
    }

    public def isConvex() : Boolean {
        return true;
    }

    public def isEmpty() : Boolean {
        return false;
    }

    public def indexOf(pt:Point) {
        if (pt.rank != 2) return -1L;
        return indexOf(pt(0), pt(1));
    }

    public def indexOf(i0:Long, i1:Long) {
        return ((i0)*(i0+1))/2 + i1;
    }

    public def min():(Int)=>Long = (i:Int)=> {
        if (i==0) return rowMin;
        else if (i==1) return colMin;
        else throw new ArrayIndexOutOfBoundsException("min: "+i+" is not a valid rank for "+this);
    };

    public def max():(Int)=>Long = (i:Int)=> {
        if (i==0) return rowMin+dim;
        else if (i==1) return colMin+dim;
        else throw new ArrayIndexOutOfBoundsException("max: "+i+" is not a valid rank for "+this);
    };

    public def size():Long { 
        return dim * (dim + 1) / 2;
    }

    public def contains(p: Point): boolean {
        if (p.rank == 2) {
            if (p(0) >= rowMin && p(0) <= (rowMin + dim)) {
                if (lower) {
                    if (p(1) >= colMin && p(1) <= p(0)) {
                        return true;
                    }
                    
                } else {
                    if (p(1) <= (colMin + dim) && p(1) >= p(0)) {
                        return true;
                    }
                }
            }
            return false;
        }
        throw new UnsupportedOperationException("contains(" + p + ")");
    }

    public def contains(r:Region(rank)): boolean {
        // TODO
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

    public def product(r: Region): Region{self!=null} {
        // TODO
        throw new UnsupportedOperationException("product()");
    }

    public def translate(v: Point(rank)): Region(rank) {
        return new TriangularRegion(rowMin + v(0), colMin + v(1), dim, lower) as Region(rank);
    }

    public def projection(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return Region.make(rowMin..(rowMin+dim));
            case 1:
                return Region.make(colMin..(colMin+dim));
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def eliminate(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return Region.make(colMin..(colMin+dim));
            case 1:
                return Region.make(rowMin..(rowMin+dim));
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def boundingBox(): Region(rank) {
        val r = Region.make(rowMin..(rowMin+dim), colMin..(colMin+dim));
        return r;
    }

    protected def computeBoundingBox(): Region(rank) {
        val r = Region.make(rowMin..(rowMin+dim), colMin..(colMin+dim));
        return r;
    }

    public def iterator(): Iterator[Point(rank)] {
        return new TriangularRegionIterator(this) as Iterator[Point(rank)];
    }

    private class TriangularRegionIterator implements Iterator[Point(2)] {
        val dim :Long;
        val lower:Boolean;
        var i:Long;
        var j:Long;

        def this(r : TriangularRegion) {
            this.dim = r.dim;
            this.lower = r.lower;
            this.i = r.rowMin;
            this.j = r.lower? r.colMin : r.colMin + r.dim;
        }

        final public def hasNext(): Boolean {
            if (i-rowMin <= dim && j-colMin <= dim) return true;
            else return false;
        }

        final public def next(): Point(2) {
            nextPoint : Point(2) = Point.make(i,j);
            if (lower) {
                if (j < (colMin + dim)) j++;
                else {
                    i++;
                    j = colMin + i;
                }
            } else {
                if (j < i) j++;
                else {
                    i++;
                    j = colMin;
                }
            }
            return nextPoint;
        } 
    }

    public def toString(): String {
        val triangleString = "triangular region " + colMin + ".." + (colMin + dim) + "," + rowMin + ".." + (rowMin + dim);
        if (lower) {
            return "lower " + triangleString;
        } else {
            return "upper " + triangleString;
        }
    }
}
