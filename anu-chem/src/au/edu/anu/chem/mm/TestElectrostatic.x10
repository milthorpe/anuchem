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

import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.util.Random;
import x10.util.ArrayList;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * A superclass for tests of electrostatic calculation methods (e.g. Direct, FMM, PME).
 * @author milthorpe
 */
public class TestElectrostatic {
    static val RANDOM_SEED = 10101110L;
    static val R = new Random(RANDOM_SEED);
    /* The maximum "noise" (random displacement) to add to particle positions from the grid. */
    static val NOISE = 0.25;

    /* side length of cubic unit cell in Angstroms */
    public val SIZE = 80.0;
    /**
     * Override this with something smaller than SIZE 
     * to produce a central cluster of particles surrounded
     * by a large empty shell.
     */
    public def sizeOfCentralCluster():Double = SIZE;
    public def boxSize():Double = SIZE;

    public def logTime(desc : String, timerIndex : Long, timer : Timer, printNewline : Boolean) {
        if (printNewline) {
            Console.OUT.printf(desc + ": %g seconds\n", (timer.mean(timerIndex) as Double) / 1e9);
        } else {
            Console.OUT.printf(desc + ": %g seconds", (timer.mean(timerIndex) as Double) / 1e9);
        }
    }

    public def logTime(desc : String, timerIndex : Long, timer : Timer) {
        logTime(desc, timerIndex, timer, true);
    }


    public def generateAtoms(numAtoms:Long) : DistArray[Rail[MMAtom]](1) {
        return generateAtoms(numAtoms, true);
    }

    /**
     * Generate an array of local Arrays of MMAtoms, one Array for each
     * place.  FMM assumes that the atoms have already been distributed.
     * Locate all particles within a small displacement from points on 
     * a cbrt(N) grid of size SIZE, centered at the origin.
     */
    public def generateAtoms(numAtoms:Long, perturb:Boolean) : DistArray[Rail[MMAtom]](1) {
        val charge = 1.0 / numAtoms; // keep charge density constant
        //Console.OUT.println("size of cluster =  " + sizeOfCentralCluster());
        val gridSize = (Math.ceil(Math.cbrt(numAtoms)) as Int);
        // assign atoms to a central cluster of size "sizeOfCentralCluster()"
        val clusterStart = -sizeOfCentralCluster() / 2.0;

        // generate particles at each place in slabs, divide space by X coordinate
        val atoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique());
        finish ateach(place in Dist.makeUnique()) {
            val chunkSize = (numAtoms / Place.MAX_PLACES) + ((numAtoms % Place.MAX_PLACES > 0) ? 1 : 0);
            val startHere = here.id * chunkSize;
            val endHere = Math.min(numAtoms, (here.id+1) * chunkSize);
            val numAtomsHere = Math.max(0L, endHere-startHere-1);
            val atomsHere = new Rail[MMAtom](numAtomsHere);
            var i:Long = 0;
            for(gridPoint in startHere..(endHere-1)) {
                val gridX = gridPoint / (gridSize * gridSize);
                val gridY = (gridPoint - (gridX * gridSize * gridSize)) / gridSize;
                val gridZ = gridPoint - (gridX * gridSize * gridSize) - (gridY * gridSize);
                val x = clusterStart + (gridX + 0.5) * (sizeOfCentralCluster() / gridSize) + (perturb ? randomNoise(): 0.0);
                val y = clusterStart + (gridY + 0.5) * (sizeOfCentralCluster() / gridSize) + (perturb ? randomNoise(): 0.0);
                val z = clusterStart + (gridZ + 0.5) * (sizeOfCentralCluster() / gridSize) + (perturb ? randomNoise(): 0.0);
                //Console.OUT.println(here.id + ": " + x + " " + y + " " + z);
                atomsHere(i++) = new MMAtom(Point3d(x, y, z), 1.0, charge);
            }
            atoms(here.id) = atomsHere;
        }

        return atoms;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    protected def getPlaceId(x : Double, y : Double, z : Double) : Int {
        return (((x + boxSize()/2.0) / boxSize()) * Place.MAX_PLACES) as Int;
    }

    /** 
     * Returns random "noise" by which to displace a particle coordinate from its
     * assigned grid point.
     */
    static def randomNoise() : Double {
        return (R.nextDouble() - 0.5) * NOISE;
    }
}
