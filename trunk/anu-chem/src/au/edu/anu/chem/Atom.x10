/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
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
    public var centre : Point3d;

    /** The symbol for this atom species. */
    public global val symbol : String;

    /** A list of atoms to which this atom is bonded. */
    private var bonds : ArrayList[Pair[BondType, Atom]]{self.at(this)};

    public def this(symbol : String, centre : Point3d) { 
        this.symbol = symbol;
        this.centre = centre;
    }

    public def this(centre : Point3d) { 
        this.symbol = "";
        this.centre = centre;
    }

    /**
     * Create a new atom with all the same values as <code>atom</code>,
     * but transferring all fields to the current place.
     */
    public def this(atom : Atom) {
        this.symbol = atom.symbol;
        this.centre = at (atom) {atom.centre};
        val bondsValRail = at(atom) {atom.bonds == null ? null : atom.bonds.toValRail()};
        if (bondsValRail != null) {
            val bonds = new ArrayList[Pair[BondType, Atom]](bondsValRail.length);
            for ((i) in 0..bondsValRail.length-1) {
                bonds.add(bondsValRail(i));
            }
            this.bonds = bonds;
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
        at (atom) { atom.addBondInternal(bondType, this); }
    }

    /**
     * Add <code>atom</code> to bond list.
     * @param bondType the type of the bond e.g. BondType.SINGLE_BOND
     * @param atom an atom to which this atom is bonded
     */
    protected def addBondInternal(bondType : BondType, atom : Atom) {
        if (bonds == null) { 
           bonds = new ArrayList[Pair[BondType, Atom]]() as ArrayList[Pair[BondType, Atom]]{self.at(this)}; 
        }
        bonds.add(Pair[BondType, Atom](bondType, atom));
    }

    /** index of this atom in a molecule */
    public var index:Int;
    public def setIndex(i:Int) { index = i; }
    public def getIndex() = index;

    public global safe def toString() : String {
        return at(this){symbol + " " + centre.i + " " + centre.j + " " + centre.k};
    }
}

