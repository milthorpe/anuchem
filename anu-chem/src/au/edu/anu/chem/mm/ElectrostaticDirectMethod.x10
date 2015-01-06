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
import x10.util.Pair;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.util.Timer;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class ElectrostaticDirectMethod implements ElectrostaticMethod {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL = 0;
    /**
     * A multi-timer for the several segments of a single invocation of
     * computePotentialAndForces, indexed by the constants above.
     */
    public val timer = new Timer(1);

    private val asyncComms = true;

    /** The charges in the simulation, divided into separate Rails for each place. */
    private val particleDataPlh:PlaceLocalHandle[ParticleData];
    private val otherAtomIndices:DistArray[Rail[Rail[Long]]](1);
    private val otherAtomData:DistArray[Rail[Rail[Double]]](1);

    /**
     * Creates a new electrostatic direct method.
     * @param particleDataPlh the particle data for each place
     */
    public def this(particleDataPlh:PlaceLocalHandle[ParticleData]) {
		this.particleDataPlh = particleDataPlh;
        this.otherAtomIndices = DistArray.make[Rail[Rail[Long]]](Dist.makeUnique(),
            (Point) => new Rail[Rail[Long]](Place.numPlaces()));
        this.otherAtomData = DistArray.make[Rail[Rail[Double]]](Dist.makeUnique(), 
            (Point) => new Rail[Rail[Double]](Place.numPlaces()));
    }
	
    public def computePotentialAndForces():Double {
        timer.start(TIMER_INDEX_TOTAL);

        val totalEnergy = finish(Reducible.SumReducer[Double]()) {
            ateach(p1 in Dist.makeUnique()) {
                offer computePotentialAndForcesLocal();
            }
        };
       
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    public def computePotentialAndForcesLocal():Double {
        val sourceId = here.id;
        var potential:Double = 0.0;
        if (asyncComms) {
            val particleData = particleDataPlh();
            val globalIndices = new Rail[Long](particleData.numAtoms(), (i:Long)=>particleData.globalIndex(i));
            val atomData = getAtomDataHere(particleData);

            // before starting computation, send my atoms to next place
            val nextPlace = Place.places().next(here);
            if (nextPlace != here) {
                @Uncounted at(nextPlace) async {
                    when(otherAtomIndices(nextPlace.id)(sourceId) == null) {
                        otherAtomIndices(nextPlace.id)(sourceId) = globalIndices;
                        otherAtomData(nextPlace.id)(sourceId) = atomData;
                    }
                }
            }

            // energy for all interactions within this place
            for (i in 0..(particleData.numAtoms()-1)) {
                val atom1Index = particleData.globalIndex(i);
                val ci = particleData.x(i);
                val xi = ci.i;
                val yi = ci.j;
                val zi = ci.k;
                val qi = particleData.atomTypes(particleData.atomTypeIndex(i)).charge;
                var fix:Double = 0.0;
                var fiy:Double = 0.0;
                var fiz:Double = 0.0;
                for (j in 0..(particleData.numAtoms()-1)) {
                    if (i==j) continue;
                    val atom2Index = particleData.globalIndex(j);
                    if (!particleData.isExcluded(atom1Index, atom2Index)) {
                        val cj = particleData.x(j);

                        val dx = cj.i-xi;
                        val dy = cj.j-yi;
                        val dz = cj.k-zi;

                        val r2 = (dx*dx + dy*dy + dz*dz);
                        val invR2 = 1.0 / r2;
                        val qq = qi * particleData.atomTypes(particleData.atomTypeIndex(j)).charge;
                        val invR = Math.sqrt(invR2);

                        val e = invR * qq;
                        val forceScaling = e * invR2;
                        potential += e;

                        val fx = forceScaling * dx;
                        val fy = forceScaling * dy;
                        val fz = forceScaling * dz;

                        fix -= fx;
                        fiy -= fy;
                        fiz -= fz;
                        //particleData.fx(j) += Vector3d(fx, fy, fz);
                    }
                }
                particleData.fx(i) = Vector3d(fix, fiy, fiz);
            }

            var target : Place = Place.places().next(nextPlace);
            var source : Place = Place.places().prev(here);
            while (source != here) {
                if (target != here) {
                    // send a set of atoms to next target place
                    val targetPlace = target;
                    @Uncounted at(targetPlace) async {
                        atomic {
                            otherAtomIndices(targetPlace.id)(sourceId) = globalIndices;
                            otherAtomData(targetPlace.id)(sourceId) = atomData;
                        }
                    }
                }
                
                var otherIndices:Rail[Long];
                var otherData:Rail[Double];
                // wait on receipt of a set of atoms from other place
                when(otherAtomIndices(here.id)(source.id) != null) {
                    otherIndices = otherAtomIndices(here.id)(source.id);
                    otherData = otherAtomData(here.id)(source.id);
                    otherAtomIndices(here.id)(source.id) = null;
                    otherAtomData(here.id)(source.id) = null;
                }

                // energy for all interactions with other atoms at other place
                for (i in 0..(particleData.numAtoms()-1)) {
                    val atom1Index = particleData.globalIndex(i);
                    val ci = particleData.x(i);
                    val xi = ci.i;
                    val yi = ci.j;
                    val zi = ci.k;
                    val qi = particleData.atomTypes(particleData.atomTypeIndex(i)).charge;
                    val fi = particleData.fx(i);
                    var fix:Double = fi.i;
                    var fiy:Double = fi.j;
                    var fiz:Double = fi.k;

                    for (j in 0..(otherIndices.size-1)) {
                        val atom2Index = otherIndices(j);
                        if (!particleData.isExcluded(atom1Index, atom2Index)) {
                            val jj = j*4;
                            val xj = otherData(jj);
                            val yj = otherData(jj+1);
                            val zj = otherData(jj+2);
                            val qj = otherData(jj+3);

                            val dx = xj-xi;
                            val dy = yj-yi;
                            val dz = zj-zi;
                            
                            val r2 = (dx*dx + dy*dy + dz*dz);
                            val invR:Double;
                            val invR2:Double;
                            invR2 = 1.0 / r2;
                            invR = Math.sqrt(invR2);
                            val qq = qi * qj;
                            val e = invR * qq;
                            potential += e;

                            val forceScaling = e * invR2;
                            fix -= forceScaling * dx;
                            fiy -= forceScaling * dy;
                            fiz -= forceScaling * dz;
                        }
                    }
                    particleData.fx(i) = Vector3d(fix, fiy, fiz);
                }
                target = Place.places().next(target);
                source = Place.places().prev(source);
            }
            
        } else {
            val particleData = particleDataPlh();
            // energy for all interactions with other atoms at other places
            finish for (p2 in Place.places()) {
                if (p2.id != sourceId) async {
                    var energyWithOther:Double = 0.0;
                    val otherPlaceAtoms = at(p2) {
                        val otherParticleData = particleDataPlh();
                        val globalIndices = new Rail[Long](otherParticleData.numAtoms(), (i:Long)=>otherParticleData.globalIndex(i));
                        val atomData = getAtomDataHere(otherParticleData);
                        new Pair[Rail[Long],Rail[Double]](globalIndices,atomData)
                    };
                    val otherIndices = otherPlaceAtoms.first;
                    val otherData = otherPlaceAtoms.second;
                    for (i in 0..(particleData.numAtoms()-1)) {
                        val atom1Index = particleData.globalIndex(i);
                        val ci = particleData.x(i);
                        val xi = ci.i;
                        val yi = ci.j;
                        val zi = ci.k;
                        val qi = particleData.atomTypes(particleData.atomTypeIndex(i)).charge;
                        val fi = particleData.fx(i);
                        var fix:Double = fi.i;
                        var fiy:Double = fi.j;
                        var fiz:Double = fi.k;

                        for (j in 0..(otherIndices.size-1)) {
                            val atom2Index = otherIndices(j);
                            if (!particleData.isExcluded(atom1Index, atom2Index)) {
                                val jj = j*4;
                                val xj = otherData(jj);
                                val yj = otherData(jj+1);
                                val zj = otherData(jj+2);
                                val qj = otherData(jj+3);

                                val dx = xj-xi;
                                val dy = yj-yi;
                                val dz = zj-zi;
                                
                                val r2 = (dx*dx + dy*dy + dz*dz);
                                val invR:Double;
                                val invR2:Double;
                                invR2 = 1.0 / r2;
                                invR = Math.sqrt(invR2);
                                val qq = qi * qj;
                                val e = invR * qq;
                                energyWithOther += e;

                                val forceScaling = e * invR2;
                                fix -= forceScaling * dx;
                                fiy -= forceScaling * dy;
                                fiz -= forceScaling * dz;
                            }
                        }
                        particleData.fx(i) = Vector3d(fix, fiy, fiz);
                    }
                    atomic potential += energyWithOther;
                }
            }

            // energy for all interactions within this place
            var energyThisPlace:Double = 0.0;
            for (i in 0..(particleData.numAtoms()-1)) {
                val atom1Index = particleData.globalIndex(i);
                for (j in 0..(i-1)) {
                    val atom2Index = particleData.globalIndex(j);
                    if (!particleData.isExcluded(atom1Index, atom2Index)) {
                        val rVec = particleData.x(j) - particleData.x(i);
                        val invR2 = 1.0 / rVec.lengthSquared();
                        val invR = Math.sqrt(invR2);
                        val atomICharge = particleData.atomTypes(particleData.atomTypeIndex(i)).charge;
                        val atomJCharge = particleData.atomTypes(particleData.atomTypeIndex(j)).charge;
                        val e = atomICharge * atomJCharge * invR;
                        energyThisPlace += 2.0 * e;
                        val pairForce = e * invR2 * rVec;
                        particleData.fx(i) -= pairForce;
		                particleData.fx(j) += pairForce;
                    }
                }
            }
            atomic potential += energyThisPlace;
        }
        return 0.5 * potential;
    }

    private def getAtomDataHere(particleData:ParticleData) {
        val atomData = new Rail[Double](particleData.numAtoms()*4);
        for (idx in 0..(particleData.numAtoms()-1)) {
            atomData(idx)   = particleData.x(idx).i;
            atomData(idx+1) = particleData.x(idx).j;
            atomData(idx+2) = particleData.x(idx).k;
            atomData(idx+3) = particleData.atomTypes(particleData.atomTypeIndex(idx)).charge;
        }
        return atomData;
    }

    public def printForces() {
        for(place in Place.places()) at(place) {
            val particleData = particleDataPlh();
            for (i in 0..(particleData.numAtoms()-1)) {
                val atomType = particleData.atomTypes(particleData.atomTypeIndex(i));
                Console.OUT.println(atomType.name + " force = " + particleData.fx(i) + " magnitude " + particleData.fx(i).length());
            }
        }
    }
}
