/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2011.
 */
package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.util.Timer;
import au.edu.anu.chem.PointCharge;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class ElectrostaticDirectMethod {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(1);

    /** The charges in the simulation, divided up into an array of ValRails, one for each place. */
    private val atoms : DistArray[Array[PointCharge](1){rect,rail}](1);

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell, divided into separate Arrays for each place
     */
    public def this(atoms : DistArray[Array[MMAtom](1){rect,rail}](1)) {
		this.atoms = DistArray.make[Array[PointCharge](1){rect,rail}](Dist.makeUnique(), 
			([i] : Point) => new Array[PointCharge](atoms(i).size(),
												(j : Int) => new PointCharge(atoms(i)(j).centre, atoms(i)(j).charge)));
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        val directEnergy = finish(SumReducer()) {
			val atoms = this.atoms; // TODO shouldn't be necessary XTENLANG-1913
            ateach ([p1] in atoms) {
                val myAtoms = atoms(p1);
                // energy for all interactions with other atoms at other places
                for ([p2] in atoms) async {
                    if (p2 != p1) { // TODO region difference
                        var energyWithOther : Double = 0.0;
                        val otherAtoms = at(atoms.dist(p2)) {atoms(p2)};
                        for ([j] in 0..(otherAtoms.size-1)) {
                            val atomJ = otherAtoms(j);
                            for ([i] in 0..(myAtoms.size-1)) {
                                val atomI = myAtoms(i);
                                energyWithOther += atomI.charge * atomJ.charge / atomJ.centre.distance(atomI.centre);
                            }
                        }
                        offer energyWithOther;
                    }
                }

                // energy for all interactions within this place
                var energyThisPlace : Double = 0.0;
                for ([i] in 0..(myAtoms.size-1)) {
					val atomI = myAtoms(i);
                    for ([j] in 0..(i-1)) {
						val atomJ = myAtoms(j);
                        energyThisPlace += 2.0 * atomI.charge * atomJ.charge / atomJ.centre.distance(atomI.centre);
                    }
                }
                offer energyThisPlace;
            }
        };
       
        timer.stop(TIMER_INDEX_TOTAL);
        return directEnergy / 2.0;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
