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
import x10.io.IOException;
import x10.util.ArrayList;
import x10.util.OptionsParser;
import x10.util.Option;

import x10x.vector.Vector3d;

import au.edu.anu.util.Timer;
import au.edu.anu.mm.uff.UniversalForceField;

/**
 * ANU-Mechanics, a molecular mechanics simulator.
 * Performs molecular mechanics using the velocity-Verlet algorithm.
 */
public class Anumm {
    /** The particle data for the simulation, of which each place holds a portion. */
    private val particleDataPlh:PlaceLocalHandle[AnummParticleData];

    /** The force field applied to the atoms in this simulation. */
    private val forceField:ForceField;

    private transient val outputFile:GromacsCoordinateFileWriter;

    private val verbose:Boolean;

    public def this(particleDataPlh:PlaceLocalHandle[AnummParticleData],
                    forceField:ForceField,
                    outputFile:GromacsCoordinateFileWriter,
                    verbose:Boolean) {
        this.particleDataPlh = particleDataPlh;
        this.forceField = forceField;
        this.outputFile = outputFile;
        this.verbose = verbose;
    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in ps
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Long) {
        Console.OUT.println("# Timestep = " + timestep + "ps, number of steps = " + numSteps);

        var step:Long = 0;
        if (verbose) printTimestep(0.0);
        val potential = forceField.getPotentialAndForces(particleDataPlh); // get initial forces

        val kinetic = getKineticEnergy();
        Console.OUT.printf("potential %10.5e kinetic %10.5e total %10.5e\n", potential, kinetic, (potential+kinetic));
        val startingEnergy = potential+kinetic;

        var energy:Double = 0.0;
        while(step < numSteps) {
            step++;
            if (verbose) printTimestep(step*timestep);
            energy = mdStep(timestep);
        }
        if (verbose) {
            outputFile.writePositions(particleDataPlh, forceField);
        }
        Console.OUT.printf("\nrelative energy drift %10.5e", (energy-startingEnergy)/energy);
    }

    private def getKineticEnergy() {
        val energy = finish(SumReducer()) {
            ateach(p in Dist.makeUnique()) {
                var kinetic:Double = 0.0;
                val particleData = particleDataPlh();
                for (i in 0..(particleData.numAtoms()-1)) {
                    val mass = particleData.atomTypes(particleData.atomTypeIndex(i)).mass;
                    kinetic += mass * particleData.dx(i).lengthSquared();
                }
                offer kinetic;
            }
        };
        return 0.5 * energy;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    private def printTimestep(time:Double) {
        Console.OUT.printf("\n%10.4f ", time);
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param dt time in ps
     * @param returns the total system energy
     */
    public def mdStep(dt:Double):Double {
        finish ateach(place in Dist.makeUnique()) {
            val particleData = particleDataPlh();

            for (i in 0..(particleData.numAtoms()-1)) {
                val invMass = 1.0 / particleData.atomTypes(particleData.atomTypeIndex(i)).mass;
                particleData.dx(i) = particleData.dx(i) + 0.5 * dt * invMass * particleData.fx(i);
                particleData.x(i) = particleData.x(i) + particleData.dx(i) * dt;
                particleData.fx(i) = Vector3d.NULL; // zero before next force calculation
            }
        }

        val potential = forceField.getPotentialAndForces(particleDataPlh);

        val kinetic = finish(SumReducer()) {
            ateach(p in Dist.makeUnique()) {
                val particleData = particleDataPlh();
                var kineticHere:Double = 0.0;
                for (i in 0..(particleData.numAtoms()-1)) {
                    val mass = particleData.atomTypes(particleData.atomTypeIndex(i)).mass;
                    val invMass = 1.0 / mass;
                    particleData.dx(i) = particleData.dx(i) + 0.5 * dt * invMass * particleData.fx(i);
                    val ke = 0.5 * mass * particleData.dx(i).lengthSquared();
                    //Console.OUT.println("kinetic for " + i + " = " + ke + " mass " + mass);
                    kineticHere += ke;
                }
                offer kineticHere;
            }
        };
        Console.OUT.printf("potential %10.5e kinetic %10.5e total %10.5e", potential, kinetic, (potential+kinetic));
        return potential+kinetic;
    }

    public static def main(args:Rail[String]) {
        val opts = new OptionsParser(args, [
            Option("h","help","this information"),
            Option("v","verbose","print out each iteration")
        ], [
            Option("c","coordinateFile","filename for structure/coordinates"),
            Option("t","timestep","timestep in picoseconds (default = 0.001)"),
            Option("n","numSteps","number of timesteps to simulate"),
            Option("p","topologyFile","GROMACS topology file"),
            Option("o","outputFile","output GROMACS coordinate file")
        ]);
        if (opts.filteredArgs().size!=0) {
            Console.ERR.println("Unexpected arguments: "+opts.filteredArgs());
            Console.ERR.println("Use -h or --help.");
            System.setExitCode(1n);
            return;
        }
        if (opts("h")) {
            Console.OUT.println(opts.usage("Usage:\n"));
            return;
        }

        val structureFileName = opts("c", null);
        if (structureFileName == null) {
            Console.OUT.println("anumm: missing coordinate file name");
            Console.OUT.println(opts.usage("Usage:\n"));
            System.setExitCode(1n);
            return;
        }
        val timestep = opts("t", 0.001);
        val numSteps = opts("n", 200n);
        val topologyFileName = opts("p", null);
        val defaultOutputFileName = structureFileName.substring(0n, structureFileName.indexOf("."))
            + "out" + numSteps + ".gro";
        val outputFileName = opts("o", defaultOutputFileName);

        var ff:ForceField = null;
        val particleDataPlh = PlaceLocalHandle.make[AnummParticleData](Place.places(), ()=> new AnummParticleData());

        val particleData = particleDataPlh();
        if (structureFileName.length() > 4) {
            val fileExtension = structureFileName.substring(structureFileName.length()-4n, structureFileName.length());
            if (fileExtension.equals(".xyz")) {
                ff = new UniversalForceField();
                try {
                    new XYZStructureFileReader(structureFileName).readParticleData(particleData, ff);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                    System.setExitCode(1n);
                    return;
                }
            } else if (fileExtension.equals(".mol") || fileExtension.equals(".sdf")) {
                ff = new UniversalForceField();
                try {
                    new MOLStructureFileReader(structureFileName).readParticleData(particleData, ff);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                    System.setExitCode(1n);
                    return;
                }
            } else if (fileExtension.equals(".gro")) {
                if (topologyFileName == null) {
                    Console.OUT.println("anumm: missing topology file\nTry anumm --help for more information.");
                    System.setExitCode(1n);
                    return;
                }
                ff = new GenericForceField(); // TODO read forcefield from topology
                try {
                    new GromacsTopologyFileReader(topologyFileName).readTopology(particleData);
                    new GromacsCoordinateFileReader(structureFileName).readParticleData(particleData);
                } catch (e:IOException) {
                    Console.ERR.println(e);
                    System.setExitCode(1n);
                    return;
                }
            }
        }

        val outputFile = new GromacsCoordinateFileWriter(outputFileName);

        Console.OUT.println("# MD for " + particleData.description + ": " + particleData.numAtoms() + " atoms");
        val anumm = new Anumm(particleDataPlh, ff, outputFile, true);
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

