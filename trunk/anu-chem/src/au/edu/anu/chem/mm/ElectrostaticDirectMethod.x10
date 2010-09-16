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
package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.util.Timer;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class ElectrostaticDirectMethod {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(6);

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private val atoms : DistArray[ValRail[MMAtom]](1);

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell
     */
    public def this(atoms : DistArray[ValRail[MMAtom]](1)) {
        this.atoms = atoms;
        // distribute all atoms to all places
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        val directEnergy = finish(SumReducer()) {
            ateach ([p1] in atoms) {
                val myAtoms = atoms(p1);
                // energy for all interactions with other atoms at other places
                foreach ([p2] in atoms) {
                    if (p2 != p1) { // TODO region difference
                        var energyWithOther : Double = 0.0;
                        val otherAtomsPacked = at(atoms.dist(p2)) {getPackedAtomsForPlace(p2)};
                        for ([j] in 0..otherAtomsPacked.length-1) {
                            for ([i] in 0..myAtoms.length-1) {
                                val myAtom = myAtoms(i);
                                energyWithOther += myAtom.charge * otherAtomsPacked(j).charge / otherAtomsPacked(j).centre.distance(myAtom.centre);
                            }
                        }
                        offer energyWithOther;
                    }
                }

                // energy for all interactions within this place
                foreach ([i] in 0..myAtoms.length-1) {
                    var energyThisPlace : Double = 0.0;
                    for ([j] in 0..i-1) {
                        energyThisPlace += 2.0 * myAtoms(i).charge * myAtoms(j).charge / myAtoms(j).centre.distance(myAtoms(i).centre);
                    }
                    offer energyThisPlace;
                }
            }
        };
       
        timer.stop(TIMER_INDEX_TOTAL);
        return directEnergy / 2.0;
    }

    static struct SumReducer implements Reducible[Double] {
        public safe def zero() = 0.0;
        public safe def apply(a:Double, b:Double) = (a + b);
    }

    /*
     * Returns all atom charges and coordinates for a place, in packed representation
     */
    private safe def getPackedAtomsForPlace(placeId : Int) : ValRail[MMAtom.PackedRepresentation] {
        val myAtoms = atoms(placeId);
        return ValRail.make[MMAtom.PackedRepresentation](myAtoms.length(), (i : Int) => myAtoms(i).getPackedRepresentation());
    }
}
