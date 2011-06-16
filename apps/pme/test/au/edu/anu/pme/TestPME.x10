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
package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
//import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed Particle Mesh Ewald implementation.
 * @author milthorpe
 */
public class TestPME extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public static def main(args : Array[String](1)) {
        var numAtoms : Int;
        var ewaldCoefficient : Double = 0.35;
        var cutoff : Double = 10.0;
        var gridSize : Int = 64;
        var splineOrder : Int = 4;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                ewaldCoefficient = Double.parseDouble(args(1));
                if (args.size > 2) {
                    cutoff = Double.parseDouble(args(2));
                    if (args.size > 3) {
                        gridSize = Int.parseInt(args(3));
                        if (args.size > 4) {
                            splineOrder = Int.parseInt(args(4));
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPME numAtoms [ewaldCoefficient] [cutoff] [gridSize] [splineOrder]");
            return;
        }

        if (splineOrder > gridSize) {
            Console.ERR.println("TestPME: splineOrder must not be greater than gridSize");
            return;
        }
        new TestPME().test(numAtoms, ewaldCoefficient, cutoff, gridSize, splineOrder);
    }

    public def test(numAtoms : Int, ewaldCoefficient : Double, cutoff : Double, gridSize : Int, splineOrder : Int) {

        val edges = [Vector3d(SIZE, 0.0, 0.0), Vector3d(0.0, SIZE, 0.0), Vector3d(0.0, 0.0, SIZE)];
        val g = gridSize;
        val gridSizes = new Array[Int](3, g);

        Console.OUT.println("Testing PME for " + numAtoms + " particles."
            //+ "\nBox edges: " + edges(0) + "," + edges(1) + "," + edges(2)
            + "\nGrid size: " + gridSize
            + " spline order: " + splineOrder 
            + " Beta: " + ewaldCoefficient 
            + " Cutoff: " + cutoff);

        val atoms = generateAtoms(numAtoms);
        val pme = new PME(edges, gridSizes, atoms, splineOrder, ewaldCoefficient, cutoff);
        pme.setup();
        val energy = pme.getEnergy();
        Console.OUT.println("energy = " + energy);

        logTime("Direct",            PME.TIMER_INDEX_DIRECT,        pme.timer);
        logTime("Self energy",       PME.TIMER_INDEX_SELF,          pme.timer);
        logTime("Grid charges",      PME.TIMER_INDEX_GRIDCHARGES,   pme.timer);
        logTime("Inverse FFT",       PME.TIMER_INDEX_INVFFT,        pme.timer);
        logTime("ThetaRecConvQ",     PME.TIMER_INDEX_THETARECCONVQ, pme.timer);
        logTime("Reciprocal energy", PME.TIMER_INDEX_RECIPROCAL,    pme.timer);
        logTime("Total",             PME.TIMER_INDEX_TOTAL,         pme.timer);
        logTime("Setup",             PME.TIMER_INDEX_SETUP,         pme.timer);
        logTime("Prefetch",          PME.TIMER_INDEX_PREFETCH,      pme.timer);
 /*
        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("cf. Direct calculation", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        // direct error comparison is only useful if there is a huge empty border around the particles
        val error = directEnergy - energy;
        Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));
*/
    }
}

