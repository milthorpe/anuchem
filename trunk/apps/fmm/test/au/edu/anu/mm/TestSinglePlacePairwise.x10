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
package au.edu.anu.mm;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests a simple pairwise electrostatic calculation at a single place.
 * @author milthorpe
 */
public class TestSinglePlacePairwise extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    private var directEnergy : Double = 0.0;

    public def this(numAtoms : Int) {
        super(numAtoms);
    }

    public static def main(args : Array[String](1)) {
        var numAtoms : Int;
        if (args.length > 0) {
            numAtoms = Int.parseInt(args(0));
        } else {
            Console.ERR.println("usage: pairwise [numAtoms]");
            return;
        }

        new TestSinglePlacePairwise(numAtoms).test();
    }

    public def test() {
        Console.OUT.println("Testing pairwise electrostatic calculation for " + numAtoms + " atoms");
        val atoms = generateAtoms();
        val myAtoms = atoms(0);

        val timer = new Timer(1);
        timer.start(0);

        finish foreach ([i] in 0..myAtoms.length-1) {
            var atomEnergy : Double = 0.0;
            for ([j] in 0..i-1) {
                atomEnergy += myAtoms(i).charge * myAtoms(j).charge / myAtoms(j).centre.distance(myAtoms(i).centre);
            }
            val atomEnergyFinal = atomEnergy;
            atomic { directEnergy += atomEnergyFinal; }
        }
        timer.stop(0);

        Console.OUT.println("directEnergy = " + directEnergy);
        logTime("Total time", 0, timer);
    }

    /** Assign all atoms to place 0. */
    safe def getPlaceId(x : Double, y : Double, z : Double) : Int {
        return 0;
    }
}

