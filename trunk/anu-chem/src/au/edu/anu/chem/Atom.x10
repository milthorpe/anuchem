/**
 * This class represents an Atom for the purpose of quantum
 * chemical and other computational chemistry codes.
 *
 * @author milthorpe
 *
 * @modified ganesh - toString() to simply print coordinate values insted of 
 *                    default (i, j , k) print.
 * @modified ganesh - moved bonding info from MMAtom here
 */
package au.edu.anu.chem;

import x10x.vector.Point3d;
import x10.util.ArrayList;
import x10.util.Pair;

public class Atom { 

    /** The location of the atomic nucleus. */
    public var centre : Point3d;

    /** The symbol for this atom species. */
    public global val symbol : String;

    /** A list of atoms to which this atom is bonded. */
    private var bonds : ArrayList[Pair[BondType, Atom]];

    public def this(symbol : String, centre : Point3d) { 
        this.symbol = symbol;
        this.centre = centre;
    }

    public def this(centre : Point3d) { 
        this.symbol = "";
        this.centre = centre;
    }

    public def getBonds() = bonds;

    /**
     * Add <code>atom</code> to bond list, and performs
     * a matching update to <code>atom</code>'s bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    public def addBond(bondType : BondType, atom : Atom) {
        addBondInternal(bondType, atom);
        at (atom) {atom.addBondInternal(bondType, this);}
    }

    /**
     * Add <code>atom</code> to bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    protected def addBondInternal(bondType : BondType, atom : Atom) {
        if (bonds == null) { bonds = new ArrayList[Pair[BondType, Atom]](); }
        bonds.add(Pair[BondType, Atom](bondType, atom));
    }

    public global safe def toString() : String {
        return at(this){symbol + " " + centre.i + " " + centre.j + " " + centre.k};
    }
}

