/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests the new distributed FMM implementation.
 * @author milthorpe
 */
public class TestFastMultipoleMethod extends TestElectrostatic {
    public def sizeOfCentralCluster() : Double = 80.0;

    public static def main(args : Array[String](1)) {
        var numAtoms:Int;
        var density:Double = 60.0;
        var dMax:Int = 5;
        var numTerms:Int = 10;
        var wellSpaced:Int = 2;
        var verbose:Boolean = false;
        var compare:Boolean = false;
        var forces:Boolean = false;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                density = Double.parseDouble(args(1));
                if (args.size > 2) {
                    dMax = Int.parseInt(args(2));
                    if (args.size > 3) {
                        numTerms = Int.parseInt(args(3));
                        if (args.size > 4) {
                            wellSpaced = Int.parseInt(args(4));
                            if (args.size > 5) {
                                if (args(5).equals("-verbose")) {
                                    verbose = true;
                                }
                                if (args.size > 6) {
                                    if (args(6).equals("-compare")) {
                                        compare = true;
                                        if (args.size > 7) {
                                            if (args(7).equals("-forces")) {
                                                forces = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: fmm numAtoms [density] [numTerms] [wellSpaced] [-verbose] [-compare] [-forces]");
            return;
        }

        new TestFastMultipoleMethod().test(numAtoms, density, dMax, numTerms, wellSpaced, verbose, compare, forces);
    }

    public def test(numAtoms:Int, density:Double, dMax:Int, numTerms:Int, wellSpaced:Int, verbose:Boolean, compare:Boolean, forces:Boolean) {
        if (verbose) {
            val numLevels = Math.max(2, (Math.log(numAtoms / density) / Math.log(8.0) + 1.0) as Int);
            Console.OUT.println("Testing FMM for " + numAtoms 
                      + " atoms, target density = " + density
                      + " numTerms = " + numTerms
                      + " wellSpaced = " + wellSpaced
                      + " numLevels = " + numLevels);

            val q = 1.0; // absolute value of charges
            val d = SIZE / Math.pow(8.0, numLevels) / 2.0;
            val e_terms:Double;
            if (wellSpaced == 1) {
                e_terms = q / ( (3.0 - Math.sqrt(3.0)) * d) * Math.pow((1.0 / Math.sqrt(3.0)), numTerms+1);
            } else if (wellSpaced == 2) {
                e_terms = q / ( (5.0 - Math.sqrt(3.0)) * d) * Math.pow((Math.sqrt(3.0) / 5.0), numTerms+1);
            } else {
                e_terms = 0.0;
            }
            Console.OUT.println("error bound in potential due to truncation: " + e_terms);
            if (wellSpaced == 1 && numTerms > 10) {
                val e_b_shift = q / ( (4.0 - 2.0*Math.sqrt(3.0)) * d) * Math.pow((Math.sqrt(3.0) / 2.0), numTerms+1);
                Console.OUT.println("error bound in potential due to B-shift: " + e_b_shift);
            }
        } else {
            Console.OUT.print(numAtoms + " atoms: ");
        }

        val atoms = generateAtoms(numAtoms);
        val fmm = new FastMultipoleMethod(density, dMax, numTerms, wellSpaced, SIZE, numAtoms);
        fmm.initialAssignment(atoms);
        val energy = fmm.calculateEnergy();
        for (i in 1..10) {
           //fmm.calculateEnergy();
        }
        
        if (verbose) {
            Console.OUT.println("energy = " + energy);

            logTime("(Tree construction)", FmmLocalData.TIMER_INDEX_TREE, fmm.localData().timer);

            logTime("Prefetch",   FmmLocalData.TIMER_INDEX_PREFETCH,  fmm.localData().timer);
            logTime("Upward",     FmmLocalData.TIMER_INDEX_UPWARD,    fmm.localData().timer);
            logTime("Downward",   FmmLocalData.TIMER_INDEX_DOWNWARD,  fmm.localData().timer);
        }

        logTime("Total",     FmmLocalData.TIMER_INDEX_TOTAL,     fmm.localData().timer, verbose);

        if (compare) {
            val direct = new ElectrostaticDirectMethod(atoms);
            val directEnergy = direct.getEnergy();
            logTime(" cf. Direct electrostatic ", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
            if (verbose) {
                //direct.printForces();
                val error = directEnergy - energy;
                Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));
                if (forces) {
                    fmm.printForces();
                }
            }
        }
    }
}

