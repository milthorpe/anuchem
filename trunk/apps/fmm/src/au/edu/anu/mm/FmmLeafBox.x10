package au.edu.anu.mm;

import x10.util.GrowableRail;

import x10x.vector.Point3d;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    public val atoms : GrowableRail[MMAtom{self.at(this)}]{self.at(this)} = new GrowableRail[MMAtom{self.at(this)}]();

    /** The U-list consists of all leaf boxes not well-separated to this box. */
    private var uList : ValRail[Point(3)];

    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : FmmBox) { 
        super(level, x, y, z, numTerms, parent);
    }

    public safe atomic def addAtom(atom : MMAtom{self.at(this)}) {
        atoms.add(atom);
    }

    public def getUList() = this.uList;

    public def setUList(uList : ValRail[Point(3)]) {
        this.uList = uList;
    }
    
    /*
     * Returns atom charges and coordinates packed into a ValRail of length 4*N
     * TODO XTENLANG-787 should return a 4*N Array
     */
    public def getPackedAtoms() : ValRail[Double] {
        if (atoms.length() > 0) {
            return ValRail.make[Double](4*atoms.length(), (i : Int) => getPackedAtomField(i));
        } else {
            return null;
        }
    }

    /**
     * Atoms are packed as [charge, centre.i, centre,j, centre.k]
     * for given i, atom index is i / 4
     */
    private def getPackedAtomField(i : Int) : Double {
        val atomIndex = i / 4;
        val atom = atoms(atomIndex);
        switch (i%4) {
            case 0:
                return atom.charge;
            case 1:
                return atom.centre.i;
            case 2:
                return atom.centre.j;
            case 3:
                return atom.centre.k;
            default:
                return Double.NaN;
        }
    }
}

