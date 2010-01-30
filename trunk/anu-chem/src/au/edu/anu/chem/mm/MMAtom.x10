package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

import au.edu.anu.chem.Atom;

/**
 * This class represents an Atom for the purpose of molecular 
 * mechanics calculations.
 * // TODO connectivity for force fields
 * @author milthorpe
 */
public class MMAtom extends Atom { 
    /** The current force acting upon this atom. */
    public var force : Vector3d;

    /** The effective charge in atomic units. */
    public global val charge : Double;

    public def this(centre : Point3d, charge : Double) {
        super(centre);
        this.charge = charge;
        this.force = Vector3d.NULL;
    } 
    
    /**
     * Create a new atom with all the same values as <code>atom</code>,
     * but transferring all fields to the current place.
     */
    public def this(atom : MMAtom) {
        super(at (atom.home) {atom.centre});
        this.charge = at (atom.home) {atom.charge};
        this.force = at (atom.home) {atom.force};
    }

    public def this(centre : Point3d) { 
        this(centre, 0.0);
    }

    public def pairEnergy(atom2 : MMAtom) : Double {
        return charge * (atom2.charge) / centre.distance(atom2.centre);
    }
}

