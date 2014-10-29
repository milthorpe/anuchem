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
package au.edu.anu.mm;

import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.util.ArrayList;
import x10.io.IOException;

import x10x.vector.Vector3d;

import au.edu.anu.util.Timer;
import au.edu.anu.mm.uff.UniversalForceField;

/**
 * ANU-Mechanics, a molecular mechanics simulator.
 * Performs molecular mechanics using the velocity-Verlet algorithm.
 */
public class Anumm {
    /** The particle data for the simulation, of which each place holds a portion. */
    private val particleDataPlh:PlaceLocalHandle[ParticleData];

    /** The force field applied to the atoms in this simulation. */
    private val forceField:ForceField;

    private val verbose:Boolean;

    public def this(particleDataPlh:PlaceLocalHandle[ParticleData],
                    forceField:ForceField,
                    verbose:Boolean) {
        this.particleDataPlh = particleDataPlh;
        this.forceField = forceField;
        this.verbose = verbose;
    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in fs (=ps/1000)
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Long) {
        Console.OUT.println("# Timestep = " + timestep + "fs, number of steps = " + numSteps);

        var step : Long = 0;
        if (verbose) {
            Console.OUT.printf("\n%.3f ", 0.0); // starting timestep
        }
        forceField.getPotentialAndForces(particleDataPlh); // get initial forces

        while(step < numSteps) {
            step++;
            if (verbose) {
                Console.OUT.printf("\n%.3f ", step*timestep);
            }
            mdStep(timestep);
        }
        if (verbose) {
            Console.OUT.println("\n# final positions:");
            finish for(place in Place.places()) at(place) {
                val particleData = particleDataPlh();
                // print final positions
                for (i in 0..(particleData.numAtoms-1)) {
                    Console.OUT.println(particleData.species(i)+ " " + particleData.x(i));
                }
            }
        }
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param timestep time in fs (=ps/1000)
     */
    public def mdStep(timestep : Double) {
        val t = timestep * 0.001;
        finish ateach(place in Dist.makeUnique()) {
            val particleData = particleDataPlh();
            for (i in 0..(particleData.numAtoms-1)) {
                val invMass = 1.0 / forceField.getAtomMass(particleData.species(i));
                particleData.dx(i) = particleData.dx(i) + 0.5 * t * invMass * particleData.fx(i);
                particleData.x(i) = particleData.x(i) + particleData.dx(i) * t;
                particleData.fx(i) = Vector3d.NULL; // zero before next force calculation
                Console.OUT.print(particleData.x(i) + " ");
            }
        }
        forceField.getPotentialAndForces(particleDataPlh);
        finish ateach(place in Dist.makeUnique()) {
            val particleData = particleDataPlh();
            for (i in 0..(particleData.numAtoms-1)) {
                val invMass = 1.0 / forceField.getAtomMass(particleData.species(i));
                particleData.dx(i) = particleData.dx(i) + 0.5 * t * invMass * particleData.fx(i);
                //Console.OUT.print(particleData.dx(i) + " ");
            }
        }
    }

    public static def main(args:Rail[String]) {
        var structureFileName : String = null;
        var timestep : Double = 0.2;
        var numSteps : Int = 200n;
        if (args.size > 0) {
            structureFileName = args(0);
            if (args.size > 1) {
                timestep = Double.parseDouble(args(1));
                if (args.size > 2) {
                    numSteps = Int.parseInt(args(2));
                }
            }
        } else {
            Console.ERR.println("usage: anumm structureFile [timestep(fs)] [numSteps]");
            return;
        }

        var ff:ForceField = null;
        val particleDataPlh = PlaceLocalHandle.make[ParticleData](Place.places(), ()=> new ParticleData());

        val particleData = particleDataPlh();
        if (structureFileName.length() > 4) {
            val fileExtension = structureFileName.substring(structureFileName.length()-4n, structureFileName.length());
            if (fileExtension.equals(".xyz")) {
                ff = new UniversalForceField();
                try {
                    new XYZStructureFileReader(structureFileName).readParticleData(particleData, ff);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                }
            } else if (fileExtension.equals(".mol") || fileExtension.equals(".sdf")) {
                ff = new UniversalForceField();
                try {
                    new MOLStructureFileReader(structureFileName).readParticleData(particleData, ff);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                }
            } else if (fileExtension.equals(".gro")) {
                ff = new GenericForceField(); // TODO read forcefield from top/itp files
                try {
                    new GromacsCoordinateFileReader(structureFileName).readParticleData(particleData, ff);
                    // TODO bonding from topology file
                    particleData.bonds = new Rail[Bond](0);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                }
            }
        }
        if (particleData == null) {
            Console.ERR.println("error: could not read structure file: " + structureFileName);
            return;
        }
        Console.OUT.println("# MD for " + particleData.description + ": " + particleData.numAtoms + " atoms");
        val anumm = new Anumm(particleDataPlh, ff, true);
        anumm.mdRun(timestep, numSteps);
        Console.OUT.println("\n# Run completed after " + numSteps + " steps");

        return;
    }

    /**
     * Partitions the molecule amongst all places, returning a distributed
     * array of Array[MMAtom], one Array for each place.  
     * MD requires that the atoms have already been distributed. 
     */
/*
    public static def assignAtoms(molecule : Molecule[MMAtom]) : DistArray[Rail[MMAtom]](1) {
        val tempAtoms = DistArray.make[ArrayList[MMAtom]](Dist.makeUnique(), (Point) => new ArrayList[MMAtom]());
        val atomList = molecule.getAtoms();
        val maxExtent = molecule.getMaxExtent();
        finish for (var i:Long = 0; i < atomList.size(); i++) {
            val atom = atomList(i);
            val p = getPlaceId(atom.centre.i, atom.centre.j, atom.centre.k, maxExtent);
            //Console.OUT.println(atom + " to " + p);
            at(Place(p)) async {
                val remoteAtom = new MMAtom(atom);
                atomic { tempAtoms(p).add(remoteAtom); }
            }
        }
        val atoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique(), ([p] : Point) => tempAtoms(p).toRail());
        return atoms;
    }
*/

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    private static def getPlaceId(x : Double, y : Double, z : Double, size : Double) : Int {
        return ((x / (size * 2) + 0.5) * Place.numPlaces()) as Int;
    }
}

