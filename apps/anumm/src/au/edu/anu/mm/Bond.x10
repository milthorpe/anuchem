/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

/**
 * A bond between two atoms.
 * @property typeIndex the index of the bond type in a list of all bond types
 */
public struct Bond(atom1Index:Long, atom2Index:Long, typeIndex:Int) { }
