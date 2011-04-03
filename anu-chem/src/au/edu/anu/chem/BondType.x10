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

    public def toString() = description;
}

