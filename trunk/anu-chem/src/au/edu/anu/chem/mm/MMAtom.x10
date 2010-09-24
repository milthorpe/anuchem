/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
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
    public val charge : Double;

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
        this.force = atom.force;
        this.velocity = atom.velocity;
    }

    public def this(centre : Point3d) { 
        this(centre, 0.0);
    }

    public def pairEnergy(atom2 : MMAtom) : Double {
        return charge * atom2.charge / centre.distance(atom2.centre);
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

    public def getPackedRepresentation() {
        return PackedRepresentation(symbol, charge, centre);
    }
}

