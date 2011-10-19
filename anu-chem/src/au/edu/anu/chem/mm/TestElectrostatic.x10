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
    public def sizeOfCentralCluster() : Double = SIZE;

    public def logTime(desc : String, timerIndex : Int, timer : Timer, printNewline : Boolean) {
        if (printNewline) {
            Console.OUT.printf(desc + ": %g seconds\n", (timer.mean(timerIndex) as Double) / 1e9);
        } else {
            Console.OUT.printf(desc + ": %g seconds", (timer.mean(timerIndex) as Double) / 1e9);
        }
    }

    public def logTime(desc : String, timerIndex : Int, timer : Timer) {
        logTime(desc, timerIndex, timer, true);
    }

    /**
     * Generate an array of local Arrays of MMAtoms, one Array for each
     * place.  FMM assumes that the atoms have already been distributed.
     * Locate all particles within a small displacement from points on 
     * a cbrt(N)-SIZE grid.
     */
    public def generateAtoms(numAtoms : Int) : DistArray[Rail[MMAtom]](1) {
        //Console.OUT.println("size of cluster =  " + sizeOfCentralCluster());
        val tempAtoms = DistArray.make[ArrayList[MMAtom]](Dist.makeUnique(), (Point) => new ArrayList[MMAtom]());
        val gridSize = (Math.ceil(Math.cbrt(numAtoms)) as Int);
        // assign atoms to a central cluster of size "sizeOfCentralCluster()"
        val clusterStart = SIZE / 2.0 - sizeOfCentralCluster() / 2.0;
        var gridPoint : Int = 0; // running total of assigned grid points
        finish for (var i : Int = 0; i < numAtoms; i++) {
            val gridX = gridPoint / (gridSize * gridSize);
            val gridY = (gridPoint - (gridX * gridSize * gridSize)) / gridSize;
            val gridZ = gridPoint - (gridX * gridSize * gridSize) - (gridY * gridSize);
            val x = clusterStart + (gridX + 0.5 + randomNoise()) * (sizeOfCentralCluster() / gridSize);
            val y = clusterStart + (gridY + 0.5 + randomNoise()) * (sizeOfCentralCluster() / gridSize);
            val z = clusterStart + (gridZ + 0.5 + randomNoise()) * (sizeOfCentralCluster() / gridSize);
            val charge = i%2==0?1:-1;
            val p = getPlaceId(x, y, z);
            async at(Place.place(p)) {
                val atom = new MMAtom(Point3d(x, y, z), charge);
                //Console.OUT.println(atom);
                atomic { tempAtoms(p).add(atom); }
            }
            gridPoint++;
        }
        val atoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique(), ([p] : Point) => (tempAtoms(p) as ArrayList[MMAtom]).toArray());
        return atoms;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    protected def getPlaceId(x : Double, y : Double, z : Double) : Int {
        return ((x / SIZE) * Place.MAX_PLACES) as Int;
    }

    /** 
     * Returns random "noise" by which to displace a particle coordinate from its
     * assigned grid point.
     */
    static def randomNoise() : Double {
        return (R.nextDouble() - 0.5) * NOISE;
    }
}
