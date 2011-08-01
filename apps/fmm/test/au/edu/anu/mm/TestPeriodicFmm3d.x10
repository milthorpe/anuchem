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

    public static def main(args : Array[String](1)) {
        var numAtoms : Int;
        var density : Double = 60.0;
        var numTerms : Int = 10;
        var numShells : Int = 10;
        var verbose : Boolean = false;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                density = Double.parseDouble(args(1));
                if (args.size > 2) {
                    numTerms = Int.parseInt(args(2));
                    if (args.size > 3) {
                        numShells = Int.parseInt(args(3));
                        if (numShells < 1) numShells = 1; // doesn't make sense to do 0 shells
                        if (args.size > 4) {
                            if (args(4).equals("-verbose")) {
                                verbose = true;
                            }
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPeriodicFmm3d numAtoms [density] [numTerms] [numShells] [-verbose]");
            return;
        }
        new TestPeriodicFmm3d().test(numAtoms, density, numTerms, numShells, verbose);
    }


    public def test(numAtoms : Int, density : Double, numTerms : Int, numShells : Int, verbose : Boolean) {
        if (verbose) {
            val numLevels = Math.max(2, (Math.log(numAtoms / density) / Math.log(8.0) + 1.0) as Int);
            Console.OUT.println("Testing Periodic FMM for " + numAtoms 
                              + " atoms, target density = " + density
                              + " numTerms = " + numTerms
                              + " numShells = " + numShells
                              + " numLevels = " + numLevels);
        } else {
            Console.OUT.print(numAtoms + " atoms: ");
        }
        

        val atoms = generateAtoms(numAtoms);
        val fmm3d = new PeriodicFmm3d(density, numTerms, Point3d(0.0, 0.0, 0.0), SIZE, numAtoms, numShells);
        fmm3d.assignAtomsToBoxes(atoms);
        val energy = fmm3d.calculateEnergy();
        
        if (verbose) {
            Console.OUT.println("energy = " + energy);

            logTime("(Tree construction)", Fmm3d.TIMER_INDEX_TREE, fmm3d.timer);

            logTime("Prefetch",   Fmm3d.TIMER_INDEX_PREFETCH,  fmm3d.timer);
            logTime("Direct",     Fmm3d.TIMER_INDEX_DIRECT,    fmm3d.timer);
            logTime("Multipole",  Fmm3d.TIMER_INDEX_MULTIPOLE, fmm3d.timer);
            logTime("Combine",    Fmm3d.TIMER_INDEX_COMBINE,   fmm3d.timer);
            logTime("Macroscopic", PeriodicFmm3d.TIMER_INDEX_MACROSCOPIC, fmm3d.timer);
            logTime("Downward",   Fmm3d.TIMER_INDEX_DOWNWARD,  fmm3d.timer);
        }

        logTime("Total",     Fmm3d.TIMER_INDEX_TOTAL,     fmm3d.timer);
        
/*
        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("cf. Direct calculation", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        val error = directEnergy - energy;
        Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));
*/
    }
}

