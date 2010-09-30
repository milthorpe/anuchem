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
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed periodic FMM 3D implementation.
 * @author milthorpe
 */
public class TestPeriodicFmm3d extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public def this(numAtoms : Int) {
        super(numAtoms);
    }

    public static def main(args : Array[String](1)) {
        var numAtoms : Int;
        var density : Double = 60.0;
        var numTerms : Int = 10;
        var numShells : Int = 10;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                density = Double.parseDouble(args(1));
                if (args.size > 2) {
                    numTerms = Int.parseInt(args(2));
                    if (args.size > 3) {
                        numShells = Int.parseInt(args(3));
                        if (numShells < 1) numShells = 1; // doesn't make sense to do 0 shells
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPeriodicFmm3d numAtoms [density] [numTerms] [numShells]");
            return;
        }
        new TestPeriodicFmm3d(numAtoms).test(density, numTerms, numShells);
    }


    public def test(density : Double, numTerms : Int, numShells : Int) {
        Console.OUT.println("Testing Periodic FMM for " + numAtoms 
                          + " atoms, target density = " + density
                          + " numTerms = " + numTerms
                          + " numShells = " + numShells);
        

        val atoms = generateAtoms();
        val fmm3d = new PeriodicFmm3d(density, numTerms, Point3d(0.0, 0.0, 0.0), SIZE, numAtoms, atoms, numShells);
        val energy = fmm3d.calculateEnergy();
        
        Console.OUT.println("energy = " + energy);

        val timer = fmm3d.timer;
        logTime("Prefetch",  Fmm3d.TIMER_INDEX_PREFETCH,  fmm3d.timer);
        logTime("Direct",    Fmm3d.TIMER_INDEX_DIRECT,    fmm3d.timer);
        logTime("Multipole", Fmm3d.TIMER_INDEX_MULTIPOLE, fmm3d.timer);
        logTime("Combine",   Fmm3d.TIMER_INDEX_COMBINE,   fmm3d.timer);
        logTime("Macroscopic", PeriodicFmm3d.TIMER_INDEX_MACROSCOPIC, fmm3d.timer);
        logTime("Transform", Fmm3d.TIMER_INDEX_TRANSFORM, fmm3d.timer);
        logTime("Far field", Fmm3d.TIMER_INDEX_FARFIELD,  fmm3d.timer);
        logTime("Total",     Fmm3d.TIMER_INDEX_TOTAL,     fmm3d.timer);
        Console.OUT.printf("Tree construction: %g seconds\n", (fmm3d.timer.total(Fmm3d.TIMER_INDEX_TREE) as Double) / 1e9);

        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("cf. Direct calculation", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        val error = directEnergy - energy;
        Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));
    }
}

