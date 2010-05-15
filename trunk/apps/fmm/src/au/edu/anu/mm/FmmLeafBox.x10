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
     * Returns atom charges and coordinates in packed representation
     */
    public def getPackedAtoms() : ValRail[MMAtom.PackedRepresentation] {
        if (atoms.length() > 0) {
            return ValRail.make[MMAtom.PackedRepresentation](atoms.length(), (i : Int) => atoms(i).getPackedRepresentation());
        } else {
            return null;
        }
    }
}

