/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */
package au.edu.anu.chem.mm;

import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests direct (N^2) electrostatic calculation.
 * @author milthorpe
 */
public class TestDirectElectrostatic extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public static def main(args : Array[String](1)) {
        var numAtoms : Int;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
        } else {
            Console.ERR.println("usage: direct numAtoms");
            return;
        }

        new TestDirectElectrostatic().test(numAtoms);
    }

    public def test(numAtoms : Int) {
        Console.OUT.println("Testing direct electrostatic for " + numAtoms + " particles.");

        val atoms = generateAtoms(numAtoms);

        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("Time", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        Console.OUT.println("energy = " + directEnergy);

    }
}

