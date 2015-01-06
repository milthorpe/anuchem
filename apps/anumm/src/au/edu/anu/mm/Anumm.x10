/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 * (C) IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.regionarray.Dist;
import x10.io.IOException;
import x10.util.ArrayList;
import x10.util.OptionsParser;
import x10.util.Option;

import x10x.vector.Vector3d;

import au.edu.anu.util.Timer;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.ParticleData;
import au.edu.anu.mm.uff.UniversalForceField;

/**
 * ANU-Mechanics, a molecular mechanics simulator.
 * Performs molecular mechanics using the velocity-Verlet algorithm.
 */
public class Anumm {
    static val ELECTRIC_CONVERSION_FACTOR = 1.389354859e2; // 1 / (4 PI e0)

    static val BOX_MARGIN = 0.00001; // to avoid particles falling outside the box

    private val particleDataPlh:PlaceLocalHandle[AnummParticleData];

    /** The force field applied to the atoms in this simulation. */
    private val forceField:ForceField;

    private transient val outputFile:GromacsCoordinateFileWriter;

    //val fmm:FastMultipoleMethod;
    val electrostatics:ElectrostaticDirectMethod;

    /* timestep length in ps */
    private val dt:Double;

    private val verbose:Boolean;

    public def this(particleDataPlh:PlaceLocalHandle[AnummParticleData],
                    forceField:ForceField,
                    outputFile:GromacsCoordinateFileWriter,
                    timestep:Double,
                    verbose:Boolean) {
        this.particleDataPlh = particleDataPlh;
        this.forceField = forceField;
        this.outputFile = outputFile;
        val particleData = particleDataPlh() as AnummParticleData;
/*
        val fmmDMax = 3n;
        val fmmNumTerms = 10n;
        val numAtoms = particleData.numAtoms();
        val fmmNumBoxes = Math.pow(8.0, fmmDMax);
        val fmmDensity = Math.ceil(numAtoms / fmmNumBoxes);
        //Console.OUT.println("fmm dMax = " + fmmDMax + " num terms = " + fmmNumTerms + " density = " + fmmDensity);
        // simulation cube centered at the origin, containing all atoms
        val boxSideLength = particleData.boxEdges.maxNorm() * (2.0 + BOX_MARGIN);
Console.OUT.println("boxSideLength " + boxSideLength);
        this.fmm = new FastMultipoleMethod(fmmDensity, fmmDMax, fmmNumTerms, 1n, boxSideLength);
*/
        val particleDataNotAnumm = PlaceLocalHandle.make[ParticleData](Place.places(), () => particleDataPlh());
        //fmm.initialAssignment(particleData.numAtoms() as Int, particleDataNotAnumm);
        electrostatics = new ElectrostaticDirectMethod(particleDataNotAnumm);
        this.dt = timestep;
        this.verbose = verbose;
    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(numSteps:Long) {
        Console.OUT.println("# Timestep = " + dt + "ps, number of steps = " + numSteps);

        var step:Long = 0;
        if (verbose) printTimestep(0.0);
        val bonded = 0.0;//forceField.computePotentialAndForces(particleDataPlh); // get initial forces
        val electrostatic = electrostatics.computePotentialAndForces() * ELECTRIC_CONVERSION_FACTOR;
        val potential = bonded + electrostatic;

        val kinetic = getKineticEnergy();
        Console.OUT.printf("potential %10.5e kinetic %10.5e total %10.5e\n", potential, kinetic, (potential+kinetic));
        val startingEnergy = potential+kinetic;

        var energy:Double = 0.0;
        while(step < numSteps) {
            step++;
            if (verbose) printTimestep(step*dt);
            energy = mdStep(step);
        }
/*
        if (verbose) {
            if (here == Place.FIRST_PLACE) {
                fmm.reduceMaxTimes();

                logTime("tree    ",  FmmLocalData.TIMER_INDEX_TREE,     FastMultipoleMethod.localData.timer);
                logTime("  (sort)   ", FmmLocalData.TIMER_INDEX_SORT,      FastMultipoleMethod.localData.timer);
                logTime("  (balance)", FmmLocalData.TIMER_INDEX_BALANCE,   FastMultipoleMethod.localData.timer);
                logTime("  (redist) ", FmmLocalData.TIMER_INDEX_REDIST,    FastMultipoleMethod.localData.timer);
                logTime("  (parents)", FmmLocalData.TIMER_INDEX_PARENTS,   FastMultipoleMethod.localData.timer);
                logTime("  (LET)    ", FmmLocalData.TIMER_INDEX_LET,       FastMultipoleMethod.localData.timer);

                logTime("Prefetch",  FmmLocalData.TIMER_INDEX_PREFETCH,  FastMultipoleMethod.localData.timer);
                logTime("Upward  ",  FmmLocalData.TIMER_INDEX_UPWARD,    FastMultipoleMethod.localData.timer);
                logTime("M2L     ",  FmmLocalData.TIMER_INDEX_M2L,       FastMultipoleMethod.localData.timer);
                logTime("Downward",  FmmLocalData.TIMER_INDEX_DOWNWARD,  FastMultipoleMethod.localData.timer);
                logTime("Total   ",  FmmLocalData.TIMER_INDEX_TOTAL,     FastMultipoleMethod.localData.timer);
            }
        }
*/
        outputFile.writePositions(particleDataPlh, forceField);
        Console.OUT.printf("\nrelative energy drift %10.5e", (energy-startingEnergy)/energy);
    }

    public def logTime(desc:String, timerIndex:Long, timer:Timer) {
        Console.OUT.printf(desc + ": %g seconds\n", (timer.mean(timerIndex) as Double) / 1e9);
    }

    private def getKineticEnergy() {
        val energy = finish(Reducible.SumReducer[Double]()) {
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

    private def printTimestep(time:Double) {
        Console.OUT.printf("\n%10.4f ", time);
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param currentStep the number of the current timestep
     * @param returns the total system energy
     */
    public def mdStep(currentStep:Long):Double {
        val electrostatic = finish(Reducible.SumReducer[Double]()) {
            ateach(p in Dist.makeUnique()) {
                //fmm.reassignAtoms(currentStep);
                val particleData = particleDataPlh();

                for (i in 0..(particleData.numAtoms()-1)) {
                    val invMass = 1.0 / particleData.atomTypes(particleData.atomTypeIndex(i)).mass;
                    particleData.dx(i) = particleData.dx(i) + 0.5 * dt * invMass * particleData.fx(i);
                    particleData.x(i) = particleData.x(i) + particleData.dx(i) * dt;
                    particleData.fx(i) = Vector3d.NULL; // zero before next force calculation
                }

                // TODO electrostatics only for non-bonded atoms
                val electrostaticLocal = electrostatics.computePotentialAndForcesLocal();

                offer electrostaticLocal * ELECTRIC_CONVERSION_FACTOR;
                for (i in 0..(particleData.numAtoms()-1)) {
                    particleData.fx(i) *= ELECTRIC_CONVERSION_FACTOR;
//Console.OUT.println("particle " + particleData.globalIndex(i) + " electrostatic force " + particleData.fx(i));
                }
            }
        };

        val bonded = 0.0;//forceField.computePotentialAndForces(particleDataPlh);

        val kinetic = finish(Reducible.SumReducer[Double]()) {
            ateach(p in Dist.makeUnique()) {
                val particleData = particleDataPlh();
                var kineticHere:Double = 0.0;
                for (i in 0..(particleData.numAtoms()-1)) {
//Console.OUT.println("particle " + particleData.globalIndex(i) + " total force " + particleData.fx(i));
                    val mass = particleData.atomTypes(particleData.atomTypeIndex(i)).mass;
                    val invMass = 1.0 / mass;
                    particleData.dx(i) = particleData.dx(i) + 0.5 * dt * invMass * particleData.fx(i);
                    val ke = 0.5 * mass * particleData.dx(i).lengthSquared();
                    kineticHere += ke;
                }
                offer kineticHere;
            }
        };
        //Console.OUT.println("electrostatic = " + electrostatic + " bonded =  " + bonded);
        val potential = electrostatic + bonded;
        Console.OUT.printf("potential %10.5e kinetic %10.5e total %10.5e\n", potential, kinetic, (potential+kinetic));
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
            particleData.generateNonBondedExclusions();
        }

        val outputFile = new GromacsCoordinateFileWriter(outputFileName);

        Console.OUT.println("# MD for " + particleData.description + ": " + particleData.numAtoms() + " atoms");
        val anumm = new Anumm(particleDataPlh, ff, outputFile, timestep, true);
        anumm.mdRun(numSteps);
        Console.OUT.println("\n# Run completed after " + numSteps + " steps");

        return;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    private static def getPlaceId(x : Double, y : Double, z : Double, size : Double) : Int {
        return ((x / (size * 2) + 0.5) * Place.numPlaces()) as Int;
    }
}

