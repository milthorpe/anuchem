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
    public static val TIMER_INDEX_TOTAL = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(1);

    private val asyncComms = true;

    /** The charges in the simulation, divided into separate Rails for each place. */
    private val atoms:PlaceLocalHandle[ParticleData];
    private val otherAtoms : DistArray[Rail[Rail[PointCharge]]](1);

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell, divided into separate Rails for each place
     */
    public def this(atoms:PlaceLocalHandle[ParticleData]) {
		this.atoms = atoms;
        this.otherAtoms = DistArray.make[Rail[Rail[PointCharge]]](Dist.makeUnique(), 
            ([p] : Point) => new Rail[Rail[PointCharge]](Place.numPlaces()));
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        val directEnergy = finish(Reducible.SumReducer[Double]()) {
            val atoms = this.atoms; // TODO shouldn't be necessary XTENLANG-1913
            val otherAtoms = this.otherAtoms; // TODO shouldn't be necessary XTENLANG-1913

            if (asyncComms) {
                ateach([p1] in Dist.makeUnique()) {
                    val myAtoms = atoms();
                    val myCharges = new Rail[PointCharge](myAtoms.numAtoms(), (i:Long)=>PointCharge(myAtoms.x(i), myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge));
                    var energyThisPlace : Double = 0.0;

                    // before starting computation, send my atoms to next place
                    val nextPlace = Place.places().next(here);
                    if (nextPlace != here) {
                        @Uncounted at(nextPlace) async {
                            atomic {
                                otherAtoms(nextPlace.id)(p1) = myCharges;
                            }
                        }
                    }

                    // energy for all interactions within this place
                    for (i in 0..(myAtoms.numAtoms()-1)) {
                        val ci = myAtoms.x(i);
                        val xi = ci.i;
                        val yi = ci.j;
                        val zi = ci.k;
                        val qi = myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge;
                        var fix:Double = 0.0;
                        var fiy:Double = 0.0;
                        var fiz:Double = 0.0;
                        for (j in 0..(myAtoms.numAtoms()-1)) {
                            if (i==j) continue;
                            val cj = myAtoms.x(j);

                            val dx = cj.i-xi;
                            val dy = cj.j-yi;
                            val dz = cj.k-zi;

                            val r2 = (dx*dx + dy*dy + dz*dz);
                            val invR2 = 1.0 / r2;
                            val qq = qi * myAtoms.atomTypes(myAtoms.atomTypeIndex(j)).charge;
                            val invR = Math.sqrt(invR2);

                            val e = invR * qq;
                            val forceScaling = e * invR2;
                            energyThisPlace += e;

                            val fx = forceScaling * dx;
                            val fy = forceScaling * dy;
                            val fz = forceScaling * dz;

                            fix += fx;
                            fiy += fy;
                            fiz += fz;
                            //myAtoms.fx(j) += Vector3d(-fx, -fy, -fz);
                        }
                        myAtoms.fx(i) = Vector3d(fix, fiy, fiz);
                    }

                    var target : Place = Place.places().next(nextPlace);
                    var source : Place = Place.places().prev(here);
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
                            for (i in 0..(myAtoms.numAtoms()-1)) {
                                val rVec = atomJ.centre - myAtoms.x(i);
                                val invR2 = 1.0 / rVec.lengthSquared();
                                val invR = Math.sqrt(invR2);
                                val atomICharge = myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge;
                                val e = atomICharge * atomJ.charge * invR;
                                energyThisPlace += e;
                                val pairForce = e * invR2 * rVec;
                                myAtoms.fx(i) += pairForce;
                            }
                        }
                        target = Place.places().next(target);
                        source = Place.places().prev(source);
                    }
                    
                    offer energyThisPlace;
                }
            } else {
                ateach([p1] in Dist.makeUnique()) {
                    val myAtoms = atoms();
                    // energy for all interactions with other atoms at other places
                    for (p2 in Place.places()) {
                        if (p2.id != p1) async {
                            var energyWithOther : Double = 0.0;
                            val otherPlaceAtoms = at(p2) {
                                val atomsHere = atoms();
                                val chargesHere = new Rail[PointCharge](atomsHere.numAtoms(), (i:Long)=>PointCharge(atomsHere.x(i), myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge));
                                chargesHere
                            };
                            for (j in 0..(otherPlaceAtoms.size-1)) {
                                val atomJ = otherPlaceAtoms(j);
                                for (i in 0..(myAtoms.numAtoms()-1)) {
                                    val rVec = atomJ.centre - myAtoms.x(i);
                                    val invR2 = 1.0 / rVec.lengthSquared();
                                    val invR = Math.sqrt(invR2);
                                    val atomICharge = myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge;
                                    val e = atomICharge * atomJ.charge * invR;
                                    energyWithOther += e;
                                    val pairForce = e * invR2 * rVec;
                                    myAtoms.fx(i) += pairForce;
                                }
                            }
                            offer energyWithOther;
                        }
                    }

                    // energy for all interactions within this place
                    var energyThisPlace : Double = 0.0;
                    for (i in 0..(myAtoms.numAtoms()-1)) {
                        for (j in 0..(i-1)) {
                            val rVec = myAtoms.x(j) - myAtoms.x(i);
                            val invR2 = 1.0 / rVec.lengthSquared();
                            val invR = Math.sqrt(invR2);
                            val atomICharge = myAtoms.atomTypes(myAtoms.atomTypeIndex(i)).charge;
                            val atomJCharge = myAtoms.atomTypes(myAtoms.atomTypeIndex(j)).charge;
                            val e = atomICharge * atomJCharge * invR;
                            energyThisPlace += 2.0 * e;
                            val pairForce = e * invR2 * rVec;
                            myAtoms.fx(i) += pairForce;
				            myAtoms.fx(j) -= pairForce;
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
        for(place in Place.places()) at(place) {
            val atomsHere = atoms();
            for (i in 0..(atomsHere.numAtoms()-1)) {
                val atomType = atomsHere.atomTypes(atomsHere.atomTypeIndex(i));
                Console.OUT.println(atomType.name + " force = " + atomsHere.fx(i) + " magnitude " + atomsHere.fx(i).length());
            }
        }
    }
}
