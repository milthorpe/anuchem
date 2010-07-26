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

