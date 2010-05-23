/*
 * This file is part of ANU Molecular Mechanics (ANUMM).
 * ANUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUMM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.chem;

/**
 * This "enum" struct type is used to differentiate among different classes of chemical bond.
 */
public struct BondType(description : String, bondOrder : Double) {
    private def this(description : String, bondOrder : Double) {
        property(description, bondOrder);
    }

    public static val WEAK_BOND = BondType("Weak bond", 0.0);
    public static val NO_BOND = BondType("No bond", 0.0);
    public static val SINGLE_BOND = BondType("Single bond", 1.0);
    public static val DOUBLE_BOND = BondType("Double bond", 2.0);
    public static val TRIPLE_BOND = BondType("Triple bond", 3.0);
    public static val QUADRUPLE_BOND = BondType("Quadruple bond", 4.0);
    public static val AROMATIC_BOND = BondType("Aromatic bond", 1.5);
    public static val AMIDE_BOND = BondType("Amide bond", 1.41);
    public static val IONIC_BOND = BondType("Ionic bond", 0.0);

    /**
     * @return true if this is a strong bond type
     */
    public def isStrongBond() : Boolean {
        return (this != BondType.NO_BOND && this != BondType.WEAK_BOND);
    }

    public global safe def toString() = description;
}

