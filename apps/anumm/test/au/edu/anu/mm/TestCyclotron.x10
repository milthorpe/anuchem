/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2013.
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.io.File;
import x10.io.FileReader;
import x10.io.EOFException;
import x10.util.ArrayList;
import x10.util.Random;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.AtomType;
import au.edu.anu.chem.mm.ParticleData;
import au.edu.anu.util.StringSplitter;

/**
 * Simulates an ion packet travelling in a circular path in a Penning trap
 * for the purposes of Fourier Transform Ion Cyclotron Resonance (FT-ICR) mass
 * spectroscopy.
 * Parameters to match simulation reported in:
 * @see Leach et al. (2010). "Comparison of Particle-In-Cell simulations with 
 *   experimentally observed frequency shifts between ions of the same mass-to-
 *   charge in Fourier Transform Ion Cyclotron Resonance Mass Spectrometry".
 *   J. Am. Soc. Mass Spectrometry 21 (2), 203-208 
 */
public class TestCyclotron {
    public static def main(args:Rail[String]) {
        if (args.size > 0) {
            val inputFile = args(0);
            var snapshotFile:String = null;
            if (args.size > 1) {
                snapshotFile = args(1);
            }
            runFromInput(inputFile, snapshotFile);
        } else {
            Console.ERR.println("usage: anum <inputFile> [snapshotFile]");
        }
    }

    /**
     * Run simulation of FT-ICR experiment according to specifications in an 
     * input file in the format described below.
     *
     * <title>
     * B <axial magnetic field in T>
     * V <quadrupolar trapping field in V>
     * edgeLength <side length of cubic Penning trap>
     * dt <timestep in ns>
     * steps <number of timesteps>
     * + fmmDMax <number of particles per lowest level box>
     * + fmmTerms <number of terms in FMM expansions>
     * * species <name> <mass> <charge> <number of ions>
     */
    public static def runFromInput(fileName:String, snapshotFileName:String) { 
        val fil = new FileReader(new File(fileName));

        var line:String = fil.readLine();

        while (line.startsWith("#")) line = fil.readLine();

        val title = line;

        var B:Double = 7.0; // magnetic field
        var V:Double = 1.0; // trapping potential
        var edgeLength:Double = 0.0508; // trap side length
        var radius:Double = 0.006; // excitation radius
        var dt:Double = 0.0;
        var steps:Int = 0n;
        var logSteps:Int = 1n;
        var fmmDensity:Double = 60.0;
        var fmmDMax:Int = 4n;
        var fmmTerms:Int = 10n;

        line = fil.readLine();
        if (line.startsWith("B")) {
            B = getDoubleParam(line);
        } else {
            throw new Exception("Invalid input: must specify magnetic field strength [B]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("V")) {
            V = getDoubleParam(line);
        } else {
            throw new Exception("Invalid input: must specify trapping potential [V]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("edgeLength")) {
            edgeLength = getDoubleParam(line);
        } else {
            throw new Exception("Invalid input: must specify trap side length [edgeLength]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("radius")) {
            radius = getDoubleParam(line);
        } else {
            throw new Exception("Invalid input: must specify excitation radius [radius]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("dt")) {
            dt = getDoubleParam(line);
        } else {
            throw new Exception("Invalid input: must specify timestep [dt]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("steps")) {
            steps = getIntParam(line);
        } else {
            throw new Exception("Invalid input: must specify number of timesteps [steps]. Next line was:\n"+line);
        }

        line = fil.readLine();
        if (line.startsWith("logSteps")) {
            logSteps = getIntParam(line);
            line = fil.readLine();
        }

        if (line.startsWith("fmmDMax")) {
            fmmDMax = getIntParam(line);
            line = fil.readLine();
        }

        if (line.startsWith("fmmTerms")) {
            fmmTerms = getIntParam(line);
            line = fil.readLine();
        }

        val particleDataPlh = PlaceLocalHandle.make[ParticleData](Place.places(), () => new ParticleData());
        val particleData = particleDataPlh();

        var totalIons:Int = 0n;
        val atomTypes = new ArrayList[AtomType]();
        val numPerSpecies = new ArrayList[Long]();
        try {
            while (line != null && line.startsWith("species")) {
                val wrd = StringSplitter.splitOnWhitespace(line);
                val name = wrd(1);
                val mass = Double.parseDouble(wrd(2));
                val charge = Int.parseInt(wrd(3));
                val numIons = Int.parseInt(wrd(4));
                totalIons += numIons;

                atomTypes.add(new AtomType(name, -1n, mass, charge));
                numPerSpecies.add(numIons);

                line = fil.readLine();
            }
        } catch (e:EOFException) {
            // no more species
        }
        if (atomTypes.isEmpty()) {
            throw new Exception("Invalid input: expected at least one species. Next line was:\n"+line);
        }
        fil.close();
        particleData.atomTypes = atomTypes.toRail();

        Console.OUT.printf("# Testing %s: cyclotron trapping potential: %2.1f V magnetic field: %6.4f T edgeLength %4.1f mm\n", title, V, B, edgeLength*1e3);
        Console.OUT.println("# species:");

        val rand = new Random(27178281L);

        var i:Int = 0n;
        for (atomTypeIndex in 0..(atomTypes.size()-1)) {
            val atomType = atomTypes(atomTypeIndex);
            val omega_c = atomType.charge * B / atomType.mass * (PenningTrap.CHARGE_MASS_FACTOR);
            val omega_z = Math.sqrt(2.0 * PenningTrap.ALPHA_PRIME * atomType.charge * V / (atomType.mass * edgeLength*edgeLength) * PenningTrap.CHARGE_MASS_FACTOR);
            val omega_plus = omega_c / 2.0 + Math.sqrt(omega_c*omega_c / 4 - omega_z*omega_z / 2);
            val nu_plus = omega_plus / (2.0 * Math.PI);

            val v = PenningTrap.CHARGE_MASS_FACTOR * B * radius * atomType.charge / atomType.mass;
            val r = atomType.mass * v / (atomType.charge * B) / PenningTrap.CHARGE_MASS_FACTOR;
            Console.OUT.printf("# %6i %10s mass %8.5f charge %i ", numPerSpecies(atomTypeIndex), atomType.name, atomType.mass, atomType.charge);
            Console.OUT.printf("omega_c = %7i rad/s omega_z = %6i rad/s omega_+ = %7i rad/s nu_+ = %7i v = %9.3g m/s\n", omega_c as Int, omega_z as Int, omega_plus as Int, nu_plus as Int, v);

            if (snapshotFileName == null) {
                for (j in 0..(numPerSpecies(atomTypeIndex)-1)) { // 
                    // distribution for each species is uniform 1mm cylinder along
                    // z dimension centred at [x=-(excitation radius), y=0, z=0]
                    val er = rand.nextDouble() * 1.0e-3;
                    val theta = rand.nextDouble() * Math.PI * 2.0;
                    val ex = Math.cos(theta) * er;
                    val ey = Math.sin(theta) * er;
                    // velocity is Maxwellian distribution with addition of velocity v in y direction
                    val ev = maxwellianVelocity(rand, atomType.mass);

                    particleData.addAtom(i++, atomTypeIndex as Int, Point3d(-radius+ex, ey, perturbation(rand, 1e-3)), Vector3d(ev.i, v+ev.j, ev.k));
                }
            }
        }

        val trap = new PenningTrap(totalIons, V, new Vector3d(0.0, 0.0, B), edgeLength, fmmDMax, fmmTerms);

        var resumeStep:Int = 0n;
        if (snapshotFileName != null) {
            resumeStep = trap.loadFromSnapshot(snapshotFileName, particleData);
        }
        trap.fmm.initialAssignment(particleDataPlh().numAtoms() as Int, particleDataPlh);

        trap.mdRun(dt, steps, resumeStep, logSteps);
    }

    private static def perturbation(rand:Random, max:Double) {
        return rand.nextDouble()*max*2.0 - max;
    }

    /** 
     * Returns a velocity vector taken from a Maxwellian distribution at 300K,
     * i.e. a vector with each component a normal distribution with mean 0 and
     * variance kT/m
     * @param rand a uniform random number generator
     * @param mass the particle mass
     */
    private static def maxwellianVelocity(rand:Random, mass:Double):Vector3d {
        val variance = 8.3144621 * 300 / mass; // (gas constant R) * 300K / (molecular mass m)
        return Vector3d(randNormal(rand, variance), randNormal(rand, variance), randNormal(rand, variance));
    }

    /**
     * Use the Box-Muller transform to generate a pseudo-random number from a normal distribution.
     * @param rand a uniform random number generator
     * @param variance distribution variance 
     * @return a random variable taken from a normal distribution with mean 0
     */
    private static def randNormal(rand:Random, variance:Double) {
        val u1 = rand.nextDouble();
        val u2 = rand.nextDouble();
        val z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        return variance * z0;
    }

    private static def getIntParam(line:String) = Int.parseInt(StringSplitter.splitOnWhitespace(line)(1));
    private static def getDoubleParam(line:String) = Double.parseDouble(StringSplitter.splitOnWhitespace(line)(1));

}

