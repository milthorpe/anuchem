/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2013.
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.chem;

import x10x.vector.Point3d;
import x10.util.ArrayList;
import x10.util.Pair;

/**
 * This class represents an Atom for the purpose of quantum
 * chemical and other computational chemistry codes.
 *
 * @author milthorpe, V. Ganesh
 */
public class Atom { 
    /** The location of the atomic nucleus. */
    public var centre:Point3d;

    /** A list of atoms to which this atom is bonded. */
    private var bonds:ArrayList[Pair[BondType, Atom]];

    /** The unique identifier for this atom species within a given simulation. */
    public val species:Int;

    /** Index of this atom within a molecule or system. */
    public var index:Long;

    public def setIndex(i:Long) { index = i; }
    public def getIndex() = index;

    public def this(species:Int, centre : Point3d) { 
        this.species = species;
        this.centre = centre;
    }

    public def this(centre : Point3d) { 
        this.species = -1n;
        this.centre = centre;
    }

    /**
     * Create a new atom with all the same values as <code>atom</code>,
     * but transferring all fields to the current place.
     */
    public def this(atom : Atom) {
        this.species = atom.species;
        this.centre = atom.centre;
        if (atom.bonds != null) {
            this.bonds = ArrayList.make(atom.bonds);
        }
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
        atom.addBondInternal(bondType, this);
    }

    /**
     * Add <code>atom</code> to bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    protected def addBondInternal(bondType : BondType, atom : Atom) {
        if (bonds == null) { 
           bonds = new ArrayList[Pair[BondType, Atom]]() as ArrayList[Pair[BondType, Atom]]; 
        }
        bonds.add(Pair[BondType, Atom](bondType, atom));
    }

    public def toString() : String {
        return String.format("%d %#10.5f %#10.5f %#10.5f", [species, centre.i, centre.j, centre.k]);
    }
}

