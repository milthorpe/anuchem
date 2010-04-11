package x10.array;

public class TriangularRegion extends BaseRegion {
    // XTENLANG-49
    static type TriangularRegion(rank:Int) = TriangularRegion{self.rank==rank};

    global val dim : Int;
    global val rowMin : Int;
    global val colMin : Int;
    global val lower : Boolean;

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

    public global def isConvex() : Boolean {
        return true;
    }

    public global def isEmpty() : Boolean {
        return false;
    }

    public global def min(): ValRail[Int] {
        return [rowMin,colMin];
    }

    public global def max(): ValRail[Int] {
        return [rowMin+dim,colMin+dim];
    }   

    public global def size() : Int { 
        return dim * (dim + 1) / 2;
    }

    public global def contains(p: Point): boolean {
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
        return new TriangularRegion(rowMin + v(0), colMin + v(1), dim, lower) as Region(rank);
    }

    public global def projection(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return rowMin..rowMin+dim;
            case 1:
                return colMin..colMin+dim;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public global def eliminate(axis: Int): Region(1) {
        switch (axis) {
            case 0:
                return colMin..colMin+dim;
            case 1:
                return rowMin..rowMin+dim;
            default:
                throw U.unsupported(this, "projection(" + axis + ")");
        }
    }

    public global def boundingBox(): Region(rank) {
        return [rowMin..rowMin+dim,colMin..colMin+dim];
    }

    protected global def computeBoundingBox(): Region(rank) {
        return [rowMin..rowMin+dim,colMin..colMin+dim];
    }

    public global def iterator(): Iterator[Point(rank)] {
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

    public global safe def toString(): String {
        val triangleString = "triangular region " + colMin + ".." + (colMin + dim) + "," + rowMin + ".." + (rowMin + dim);
        if (lower) {
            return "lower " + triangleString;
        } else {
            return "upper " + triangleString;
        }
    }
}
