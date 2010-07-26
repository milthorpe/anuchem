/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.Atom;
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

    public def this(symbol : String, centre : Point3d, charge : Double) {
        super(symbol, centre);
        this.charge = charge;
        this.force = Vector3d.NULL;
        this.velocity = Vector3d.NULL;
    }

    public def this(centre : Point3d, charge : Double) {
        super(centre);
        this.charge = charge;
        this.force = Vector3d.NULL;
        this.velocity = Vector3d.NULL;
    } 
    
    /**
     * Create a new atom with all the same values as <code>atom</code>,
     * but transferring all fields to the current place.
     */
    public def this(atom : MMAtom) {
        super(atom);
        this.charge = atom.charge;
        this.force = at (atom) {atom.force};
        this.velocity = at (atom) {atom.velocity};
    }

    public def this(centre : Point3d) { 
        this(centre, 0.0);
    }

    public def pairEnergy(atom2 : MMAtom) : Double {
        return charge * atom2.charge / centre.distance(at(atom2){atom2.centre});
    }

    /**
     * A packed representation of an MMAtom, for use in
     * transferring bulk Atom data between places.
     * Does not include force or velocity, as these
     * are not commonly required at remote places.
     */
    public static struct PackedRepresentation(symbol : String, charge : Double, centre : Point3d) { 
        public def this(symbol : String, charge : Double, centre : Point3d) {
            property(symbol, charge, centre);
        }
    }

    public safe def getPackedRepresentation() {
        return PackedRepresentation(symbol, charge, centre);
    }
}

