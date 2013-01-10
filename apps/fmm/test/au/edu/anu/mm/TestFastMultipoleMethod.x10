/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012-2013.
 */
package au.edu.anu.mm;

import x10.compiler.Ifdef;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;
//import edu.utk.cs.papi.PAPI;

/**
 * Tests the new distributed FMM implementation.
 * @author milthorpe
 */
public class TestFastMultipoleMethod extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public static def main(args : Array[String](1)) {
        var numAtoms:Int;
        var dMax:Int = 3;
        var numTerms:Int = 12;
        var wellSpaced:Int = 1;
        var verbose:Boolean = false;
        var compare:Boolean = false;
        var rms:Boolean = false;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                dMax = Int.parseInt(args(1));
                if (args.size > 2) {
                    numTerms = Int.parseInt(args(2));
                    if (args.size > 3) {
                        wellSpaced = Int.parseInt(args(3));
                        if (args.size > 4) {
                            if (args(4).equals("-verbose")) {
                                verbose = true;
                            }
                            if (args.size > 5) {
                                if (args(5).equals("-compare")) {
                                    compare = true;
                                    if (args.size > 6) {
                                        if (args(6).equals("-rms")) {
                                            rms = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: fmm numAtoms [dMax] [numTerms] [wellSpaced] [-verbose] [-compare] [-rms]");
            return;
        }

        new TestFastMultipoleMethod().test(numAtoms, dMax, numTerms, wellSpaced, verbose, compare, rms);
    }

    public def test(numAtoms:Int, dMax:Int, numTerms:Int, wellSpaced:Int, verbose:Boolean, compare:Boolean, rms:Boolean) {
        val numBoxes = Math.pow(8.0, dMax);
        val density = Math.ceil(numAtoms / numBoxes);
        if (verbose) {
            Console.OUT.println("Testing FMM for " + numAtoms 
                      + " atoms, target density = " + density
                      + " dMax = " + dMax
                      + " numTerms = " + numTerms
                      + " wellSpaced = " + wellSpaced);

            val q = 1.0; // absolute value of charges
            val d = SIZE / Math.pow(2.0, dMax) / 2.0;
            val e_terms:Double;
            if (wellSpaced == 1) {
                e_terms = q / ( (3.0 - Math.sqrt(3.0)) * d) * Math.pow((1.0 / Math.sqrt(3.0)), numTerms+1);
            } else if (wellSpaced == 2) {
                e_terms = q / ( (5.0 - Math.sqrt(3.0)) * d) * Math.pow((Math.sqrt(3.0) / 5.0), numTerms+1);
            } else {
                e_terms = 0.0;
            }
            Console.OUT.println("max error in potential due to truncation: " + e_terms);
            if (wellSpaced == 1 && numTerms > 10) {
                val e_b_shift = q / ( (4.0 - 2.0*Math.sqrt(3.0)) * d) * Math.pow((Math.sqrt(3.0) / 2.0), numTerms+1);
                Console.OUT.println("max error in potential due to B-shift: " + e_b_shift);
            }
        } else {
            Console.OUT.printf("N: %8i ", numAtoms);
        }

        val atoms = generateAtoms(numAtoms, false);
        val fmm = new FastMultipoleMethod(density, dMax, numTerms, wellSpaced, SIZE, verbose);
        fmm.initialAssignment(numAtoms, atoms);

        finish ateach(place in Dist.makeUnique()) {
            for (timerId in 4..9) {
                FastMultipoleMethod.localData.timer.clear(timerId);
            }
          fmm.reassignAtoms(0);
        }
        //fmm.countOctants();

/* 
        val papi = new PAPI();
@Ifdef("__PAPI__")
{
        papi.initialize();
        papi.countFlops();
        papi.start();
}
*/
        val energy = fmm.calculateEnergy();


        for (i in 1..9) {
            finish ateach(place in atoms) {
                fmm.reassignAtoms(i);
                fmm.calculateEnergyLocal();
            }
        }

/*
@Ifdef("__PAPI__")
{
        papi.stop();
        papi.printFlops();
        papi.shutDown();
}
*/

        fmm.reduceMaxTimes();
        
        if (verbose) {
            Console.OUT.println("energy = " + energy);

            logTime("tree    ",  FmmLocalData.TIMER_INDEX_TREE,     FastMultipoleMethod.localData.timer);
            logTime("  (sort)   ", FmmLocalData.TIMER_INDEX_SORT,      FastMultipoleMethod.localData.timer);
            logTime("  (balance)", FmmLocalData.TIMER_INDEX_BALANCE,   FastMultipoleMethod.localData.timer);
            logTime("  (redist) ", FmmLocalData.TIMER_INDEX_REDIST,    FastMultipoleMethod.localData.timer);
            logTime("  (parents)", FmmLocalData.TIMER_INDEX_PARENTS,   FastMultipoleMethod.localData.timer);
            logTime("  (LET)    ", FmmLocalData.TIMER_INDEX_LET,       FastMultipoleMethod.localData.timer);

            logTime("Prefetch",  FmmLocalData.TIMER_INDEX_PREFETCH,  FastMultipoleMethod.localData.timer);
            logTime("Upward  ",  FmmLocalData.TIMER_INDEX_UPWARD,    FastMultipoleMethod.localData.timer);
            logTime("M2L     ",  FmmLocalData.TIMER_INDEX_M2L,       FastMultipoleMethod.localData.timer);
            logTime("Downward",  FmmLocalData.TIMER_INDEX_DOWNWARD,  FastMultipoleMethod.localData.timer);
            logTime("Total   ",  FmmLocalData.TIMER_INDEX_TOTAL,     FastMultipoleMethod.localData.timer);
        } else {
            Console.OUT.printf("p: %2i D_max: %2i time (s): %9.5f tree (s): %9.5f", numTerms, dMax, (FastMultipoleMethod.localData.timer.mean(FmmLocalData.TIMER_INDEX_TOTAL) as Double) / 1e9, (FastMultipoleMethod.localData.timer.mean(FmmLocalData.TIMER_INDEX_TREE) as Double) / 1e9);
        }

        if (compare) {
            if (rms) {
                fmm.printRMSErrors();
            }
            val direct = new ElectrostaticDirectMethod(atoms);
            val directEnergy = direct.getEnergy();
            logTime(" cf. Direct electrostatic ", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
            val error = Math.abs(directEnergy - energy) / Math.abs(energy);
            if (verbose) {
                //direct.printForces();

                Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + error);
            } else {
                Console.OUT.printf("E err: %.2G", error);
            }
        }
        Console.OUT.println();
    }
}

