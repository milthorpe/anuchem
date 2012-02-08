package x10.array;

public class TriangularRegion extends Region {
    // XTENLANG-49
    static type TriangularRegion(rank:Int) = TriangularRegion{self.rank==rank};

    val dim : Int;
    val rowMin : Int;
    val colMin : Int;
    val lower : Boolean;

    /**
     * Constructs a new TriangularRegion 
     * @param p the dimension of the triangular array
     */
    public def this(rowMin: Int, colMin: Int, size: Int, lower: Boolean): TriangularRegion(2) {
        super(2, false, true);
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
        if (pt.rank != 2) return -1;
        return ((pt(0) * pt(0)) / 2) + pt(1);
    }

    public def min():(int)=>int = (i:int)=> {
        if (i==0) return rowMin;
        else if (i==1) return colMin;
        else throw new ArrayIndexOutOfBoundsException("min: "+i+" is not a valid rank for "+this);
    };

    public def max():(int)=>int = (i:int)=> {
        if (i==0) return rowMin+dim;
        else if (i==1) return colMin+dim;
        else throw new ArrayIndexOutOfBoundsException("max: "+i+" is not a valid rank for "+this);
    };

    public def size() : Int { 
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
                return rowMin..(rowMin+dim);
            case 1:
                return colMin..(colMin+dim);
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def eliminate(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return colMin..(colMin+dim);
            case 1:
                return rowMin..(rowMin+dim);
            default:
                throw new UnsupportedOperationException("projection(" + axis + ")");
        }
    }

    public def boundingBox(): Region(rank) {
        val r = (rowMin..(rowMin+dim) * colMin..(colMin+dim)) as Region(2);
        return r;
    }

    protected def computeBoundingBox(): Region(rank) {
        val r = (rowMin..(rowMin+dim) * colMin..(colMin+dim)) as Region(2);
        return r;
    }

    public def iterator(): Iterator[Point(rank)] {
        return new TriangularRegionIterator(this) as Iterator[Point(rank)];
    }

    private class TriangularRegionIterator implements Iterator[Point(2)] {
        val dim : Int;
        val lower : Boolean;
        var i : Int;
        var j : Int;

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
