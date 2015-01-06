/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011-2013.
 */
package au.edu.anu.chem.mm;

import x10.compiler.Ifdef;

import au.edu.anu.util.Timer;

import edu.utk.cs.papi.PAPI;

/**
 * Tests direct (N^2) electrostatic calculation.
 * @author milthorpe
 */
public class TestDirectElectrostatic extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    transient var papi:PAPI=new PAPI(); // @Ifdef("__PAPI__") // XTENLANG-3132

    public static def main(args:Rail[String]) {
        var numAtoms:Long;
        if (args.size > 0) {
            numAtoms = Long.parseLong(args(0));
        } else {
            numAtoms = 10000;
        }

        new TestDirectElectrostatic().test(numAtoms);
    }

    public def test(numAtoms:Long) {
        Console.OUT.println("Testing direct electrostatic for " + numAtoms + " particles.");

        val atoms = generateAtoms(numAtoms);

        val direct = new ElectrostaticDirectMethod(atoms);

        @Ifdef("__PAPI__") {
        papi.initialize();
        papi.countFlops();
        papi.start();
        }

        val directEnergy = direct.computePotentialAndForces();

        @Ifdef("__PAPI__"){
        papi.stop();
        papi.printFlops();
        papi.shutDown();
        }

        logTime("Time", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        Console.OUT.println("energy = " + directEnergy);
    }
}

