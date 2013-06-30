/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.chem.mm;

import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.compiler.Uncounted;
import x10.util.Team;

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

    private val asyncComms = true;

    /** The charges in the simulation, divided into separate Rails for each place. */
    private val atoms : DistArray[Rail[MMAtom]](1);
    private val otherAtoms : DistArray[Rail[Rail[PointCharge]]](1);

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell, divided into separate Rails for each place
     */
    public def this(atoms : DistArray[Rail[MMAtom]](1)) {
		this.atoms = atoms;
        this.otherAtoms = DistArray.make[Rail[Rail[PointCharge]]](Dist.makeUnique(), 
            ([p] : Point) => new Rail[Rail[PointCharge]](Place.MAX_PLACES));
    }

    public def expectationValue(twoParticleFunction:(a:MMAtom,b:MMAtom) => Double):Double {
        var total:Double = 0.0;
        val a = atoms(here.id);
        for (i in 0..(a.size-1)) {
            for (j in 0..(a.size-1)) {
                if (i != j) {
                    total += twoParticleFunction(a(i), a(j));
                }
            }
        }
        val N = a.size;
        return total / (N*(N-1));
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

//        val radialDistribution = (a:PointCharge,b:PointCharge) => b.centre.distance(a.centre);

        val directEnergy = finish(SumReducer()) {
            val atoms = this.atoms; // TODO shouldn't be necessary XTENLANG-1913
            val otherAtoms = this.otherAtoms; // TODO shouldn't be necessary XTENLANG-1913

            if (asyncComms) {
                ateach([p1] in atoms) {
                    val myAtoms = atoms(p1);
                    val myCharges = new Rail[PointCharge](myAtoms.size, (i:Long)=>PointCharge(myAtoms(i).centre, myAtoms(i).charge));
                    var energyThisPlace : Double = 0.0;

                    // before starting computation, send my atoms to next place
                    val nextPlace = here.next();
                    if (nextPlace != here) {
                        @Uncounted at(nextPlace) async {
                            atomic {
                                otherAtoms(nextPlace.id)(p1) = myCharges;
                            }
                        }
                    }

                    // energy for all interactions within this place
                    for (i in 0..(myAtoms.size-1)) {
			            val atomI = myAtoms(i);
                        for (j in 0..(i-1)) {
				            val atomJ = myAtoms(j);
                            val rVec = atomJ.centre - atomI.centre;
                            val invR2 = 1.0 / rVec.lengthSquared();
                            val invR = Math.sqrt(invR2);
                            val e = atomI.charge * atomJ.charge * invR;
                            energyThisPlace += 2.0 * e;
                            val pairForce = e * invR2 * rVec;
                            atomI.force += pairForce;
				            atomJ.force -= pairForce;
                        }
                    }

                    var target : Place = nextPlace.next();
                    var source : Place = here.prev();
                    while (source != here) {
                        if (target != here) {
                            // send a set of atoms to next target place
                            val targetPlace = target;
                            @Uncounted at(targetPlace) async {
                                atomic {
                                    otherAtoms(targetPlace.id)(p1) = myCharges;
                                }
                            }
                        }
                        
                        // wait on receipt of a set of atoms from other place
                        when(otherAtoms(here.id)(source.id) != null);

                        // energy for all interactions with other atoms at other place
                        val other = otherAtoms(here.id)(source.id);
                        for (j in 0..(other.size-1)) {
                            val atomJ = other(j);
                            for (i in 0..(myAtoms.size-1)) {
                                val atomI = myAtoms(i);
                                val rVec = atomJ.centre - atomI.centre;
                                val invR2 = 1.0 / rVec.lengthSquared();
                                val invR = Math.sqrt(invR2);
                                val e = atomI.charge * atomJ.charge * invR;
                                energyThisPlace += e;
                                val pairForce = e * invR2 * rVec;
                                atomI.force += pairForce;
                            }
                        }
                        target = target.next();
                        source = source.prev();
                    }
                    
                    offer energyThisPlace;
                }
            } else {
                ateach([p1] in atoms) {
                    val myAtoms = atoms(p1);
                    // energy for all interactions with other atoms at other places
                    for ([p2] in atoms) async {
                        if (p2 != p1) { // TODO region difference
                            var energyWithOther : Double = 0.0;
                            val otherPlaceAtoms = at(atoms.dist(p2)) {
                                val atomsHere = atoms(p2);
                                val chargesHere = new Rail[PointCharge](atomsHere.size, (i:Long)=>PointCharge(atomsHere(i).centre, atomsHere(i).charge));
                                chargesHere
                            };
                            for (j in 0..(otherPlaceAtoms.size-1)) {
                                val atomJ = otherPlaceAtoms(j);
                                for (i in 0..(myAtoms.size-1)) {
                                    val atomI = myAtoms(i);
                                    val rVec = atomJ.centre - atomI.centre;
                                    val invR2 = 1.0 / rVec.lengthSquared();
                                    val invR = Math.sqrt(invR2);
                                    val e = atomI.charge * atomJ.charge * invR;
                                    energyWithOther += e;
                                    val pairForce = e * invR2 * rVec;
                                    atomI.force += pairForce;
                                }
                            }
                            offer energyWithOther;
                        }
                    }

                    // energy for all interactions within this place
                    var energyThisPlace : Double = 0.0;
                    for (i in 0..(myAtoms.size-1)) {
					    val atomI = myAtoms(i);
                        for (j in 0..(i-1)) {
						    val atomJ = myAtoms(j);
                            val rVec = atomJ.centre - atomI.centre;
                            val invR2 = 1.0 / rVec.lengthSquared();
                            val invR = Math.sqrt(invR2);
                            val e = atomI.charge * atomJ.charge * invR;
                            energyThisPlace += 2.0 * e;
                            val pairForce = e * invR2 * rVec;
                            atomI.force += pairForce;
				            atomJ.force -= pairForce;
                        }
                    }
                    offer energyThisPlace;
                }
            }
        };
       
        timer.stop(TIMER_INDEX_TOTAL);
        return 0.5 * directEnergy;
    }

    public def printForces() {
        for(p1 in atoms) at(atoms.dist(p1)) {
            val atomsHere = atoms(p1);
            for (i in 0..(atomsHere.size-1)) {
                val atom = atomsHere(i);
                Console.OUT.println(atom.species + " force = " + atom.force + " magnitude " + atom.force.length());
            }
        }
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
