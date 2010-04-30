package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.Atom;
import x10.util.ArrayList;
import x10.util.Pair;

/**
 * This class represents an Atom for the purpose of molecular 
 * mechanics calculations.
 * // TODO connectivity for force fields
 * @author milthorpe
 */
public class MMAtom extends Atom { 
    /** The current force acting upon this atom. */
    public var force : Vector3d;
    
    /** The current velocity of this atom. */
    public var velocity : Vector3d;

    /** The effective charge in atomic units. */
    public global val charge : Double;

    /** The mass of this atom in atomic units. */
    public global val mass : Double;

    /** A list of atoms to which this atom is bonded. */
    private var bonds : ArrayList[Pair[BondType, MMAtom]];

    public def this(symbol : String, centre : Point3d, mass : Double, charge : Double) {
        super(symbol, centre);
        this.mass = mass;
        this.charge = charge;
        this.force = Vector3d.NULL;
        this.velocity = Vector3d.NULL;
    }

    public def this(centre : Point3d, mass : Double, charge : Double) {
        super(centre);
        this.mass = mass;
        this.charge = charge;
        this.force = Vector3d.NULL;
        this.velocity = Vector3d.NULL;
    } 
    
    /**
     * Create a new atom with all the same values as <code>atom</code>,
     * but transferring all fields to the current place.
     */
    public def this(atom : MMAtom) {
        super(at (atom) {atom.centre});
        this.mass = atom.mass;
        this.charge = atom.charge;
        this.force = at (atom) {atom.force};
        this.velocity = at (atom) {atom.velocity};
    }

    public def this(centre : Point3d) { 
        this(centre, 0.0, 0.0);
    }

    public def pairEnergy(atom2 : MMAtom) : Double {
        return charge * atom2.charge / centre.distance(at(atom2){atom2.centre});
    }

    public def getBonds() = bonds;

    /**
     * Add <code>atom</code> to bond list, and performs
     * a matching update to <code>atom</code>'s bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    public def addBond(bondType : BondType, atom : MMAtom) {
        addBondInternal(bondType, atom);
        at (atom) {atom.addBondInternal(bondType, this);}
    }

    /**
     * Add <code>atom</code> to bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    def addBondInternal(bondType : BondType, atom : MMAtom) {
        if (bonds == null) { bonds = new ArrayList[Pair[BondType, MMAtom]](); }
        bonds.add(Pair[BondType, MMAtom](bondType, atom));
    }
}

