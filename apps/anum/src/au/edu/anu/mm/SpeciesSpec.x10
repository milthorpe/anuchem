/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2013.
 */
package au.edu.anu.mm;

public class SpeciesSpec {
    public val name:String;
    public val mass:Double;
    public val charge:Double;
    public val number:Int;
    public def this(name:String, mass:Double, charge:Double, number:Int) {
        this.name = name;
        this.mass = mass;
        this.charge = charge;
        this.number = number;
    }
}
