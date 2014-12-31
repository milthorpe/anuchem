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
 * This struct represents the bond stretch parameters for a pair of atom types.
 * @species the species of the first atom
 * @species the species of the second atom
 * @property idealRadius the equilibrium bond length
 * @property forceConstant the bond stretch force constant
 */
public struct BondStretchParameters(species1:Int, species2:Int, idealRadius:Double, forceConstant:Double) { }
