package x10.array;

public value class TriangularRegion extends BaseRegion{rank==2} {
    val dim : Int;
    val rowMin : Int;
    val colMin : Int;
    val lower : Boolean;

    /**
     * Constructs a new TriangularRegion 
     * @param p the dimension of the triangular array
     */
    public def this(rowMin: Int, colMin: Int, size: Int, lower: Boolean) {
        super(size, false, true);
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

    public def min(): ValRail[Int] {
        return [rowMin,colMin];
    }

    public def max(): ValRail[Int] {
        return [rowMin+dim,colMin+dim];
    }   

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

    public def projection(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return rowMin..rowMin+dim;
            case 1:
                return colMin..colMin+dim;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public def eliminate(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return colMin..colMin+dim;
            case 1:
                return rowMin..rowMin+dim;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public def boundingBox(): Region(rank) {
        return [rowMin..rowMin+dim,colMin..colMin+dim];
    }

    public def computeBoundingBox(): Region(rank) {
        return [rowMin..rowMin+dim,colMin..colMin+dim];
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
            nextPoint : Point(2) = [i,j];
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

}
