/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.io.File;
import x10.io.FileWriter;
import x10.io.IOException;
import x10.io.Printer;
import x10.util.ArrayList;
import x10.util.Team;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.StatisticalTimer;
import edu.mit.fftw.FFTW;

/**
 * A simulation of charged ions in cyclotron motion in a cubic Penning trap.
 * Equations of motion are integrated using the Boris scheme.
 * Lengths are in nm, hence the constant factors 1e-9.
 * @see Shenheng Guan and Alan G Marshall (1995)
 *   "Ion traps for Fourier transform ion cyclotron resonance mass spectrometry:
 *   principles and design of geometric and electric configurations"
 *   Int. J. Mass Spectrometry and Ion Processes 146/147 261-296
 * @see Seung-Jin Han and Seung Koo Shin (1997)
 *   "Space-Charge Effects and Fourier transform ion cyclotron resonance signals:
 *   experimental observations and three-dimensional trajectory simulations"
 *   J. Am. Soc. Mass Spectrometry 8 (4), 319-326 
 * @see http://physics.nist.gov/cuu/Constants/index.html for physical constants
 */
public class PenningTrap {
    public static val LENGTH_FACTOR = 1.0e-9;
    static val CHARGE_FACTOR = 1.602176565e-19; // conversion from e to C
    public static val CHARGE_MASS_FACTOR = 9.64853365e7; // conversion of q/m from e/Da to C/kg
    public static val MASS_CHARGE_FACTOR = 1.036426919e-8; // conversion of m/q from Da/e to kg/C
    static val ALPHA_PRIME = 2.77373; // geometric factor for a cubic trap (Guan and Marshall eq. 59)
    static val BETA_PRIME = 0.72167; // electric field constant for detection/excition (Guan and Marshall eq. 66)
    static val COULOMB_FACTOR = 1.439964485e-9; // CHARGE_FACTOR / (4 PI e0)
    private var numAtoms:Int;

    /** Side length of cubic trap. */
    val edgeLength:Double;

    /** Normalization factor for electric field. */
    val eNorm:Double;

    /** Electric field factor (BETA_PRIME / edgeLength) for detection and excitation. */
    val eFieldNorm:Double;

    val fmm:FastMultipoleMethod;

    private val V_T:Double; // trapping potential, in V

    /** The static homogeneous magnetic field B, in Teslas. */
    public val B:Vector3d;
    /** The scalar magnitude of B. */
    private val magB:Double;

    /** The maximum calculated error in particle x displacement. */
    private var maxErrorX:Double = 0.0;

    /** 
     * Creates a new Penning trap containing the given atoms.
     * @param trappingPotential the axial confinement potential in V applied to the end plates
     * @param magneticField the radial confinement field in T
     * @param edgeLength the side length of the cubic trap
     */
    public def this(numAtoms:Int,
                    atoms:DistArray[Rail[MMAtom]](1),
                    trappingPotential:Double,
                    magneticField:Vector3d,
                    edgeLength:Double,
                    fmmDMax:Int,
                    fmmNumTerms:Int) {
        this.numAtoms = numAtoms;
        val wellSeparatedParam = 2;
        val fmmNumBoxes = Math.pow(8.0, fmmDMax);
        val fmmDensity = Math.ceil(numAtoms / fmmNumBoxes);
        Console.OUT.println("fmm num levels = " + fmmDMax + " density = " + fmmDensity);
        this.fmm = new FastMultipoleMethod(fmmDensity, fmmDMax, fmmNumTerms, wellSeparatedParam, edgeLength);
        fmm.initialAssignment(numAtoms, atoms);
        this.V_T = trappingPotential;
        this.B = magneticField;
        this.edgeLength = edgeLength;
        this.eNorm = V_T * PenningTrap.ALPHA_PRIME / (edgeLength * edgeLength);
        this.eFieldNorm = -BETA_PRIME / edgeLength;
        this.magB = magneticField.magnitude();

    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in ns
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Int, logSteps:Int) {
        Console.OUT.println("# Timestep = " + timestep + "ns, number of steps = " + numSteps);
        val timer = new StatisticalTimer(3);
 
        SystemProperties.printHeader();

        // initialise position files
        val startPosHeaderPrinter = new Printer(new FileWriter(new File("positions_0.dat"), false));
        startPosHeaderPrinter.println("# positions at time 0");

        val endTime = (timestep * numSteps) as Int;
        val endPosHeaderPrinter = new Printer(new FileWriter(new File("positions_" + endTime + ".dat"), false));
        endPosHeaderPrinter.println("# positions at time " + endTime);

        // intitalise energy file
        val energiesHeaderPrinter = new Printer(new FileWriter(new File("energies_" + endTime + ".dat"), false));
        endPosHeaderPrinter.println("# energies at time " + endTime);

        finish ateach(place in Dist.makeUnique()) {
            var step:Int = 0;
            val props = new SystemProperties();
            printStartPositions(props);

            val current = new Array[Double](numSteps);

            val timeStepSecs = timestep * 1.0e-9;
            while(step < numSteps) {
                step++;
                val accumProps = (step % logSteps == 0);
                mdStepLocal(step, timeStepSecs, current, accumProps, props, timer);
                if (accumProps) {
                    reduceAndPrintProperties(timestep * step, props);
                } else {
                    Team.WORLD.barrier(here.id); // TODO is this really required?
                }
            }

            // gather current data
            Team.WORLD.allreduce[Double](here.id, current, 0, current, 0, current.size, Team.ADD);
            if (here == Place.FIRST_PLACE) {
                printCurrent(timestep, current);
                printMassSpectrum(timestep, current);
            }

            printEndPositions(timestep * step);
            printEnergies(timestep * step);

            Team.WORLD.allreduce[Long](here.id, timer.count, 0, timer.count, 0, timer.count.size, Team.ADD);
            Team.WORLD.allreduce[Long](here.id, timer.total, 0, timer.total, 0, timer.total.size, Team.ADD);
            Team.WORLD.allreduce[Long](here.id, timer.min, 0, timer.min, 0, timer.min.size, Team.MIN);
            Team.WORLD.allreduce[Long](here.id, timer.max, 0, timer.max, 0, timer.max.size, Team.MAX);
            Team.WORLD.allreduce[Double](here.id, timer.sumOfSquares, 0, timer.sumOfSquares, 0, timer.sumOfSquares.size, Team.ADD);
            if (here == Place.FIRST_PLACE) {
                timer.printSeconds();
                Console.OUT.println("fmm:");
                timer.printSeconds(1);
                Console.OUT.println("reassign:");
                timer.printSeconds(2);
            }
        }

    }

    private def printStartPositions(props:SystemProperties) {
        val startPosPrinter = new Printer(new FileWriter(new File("positions_0.dat"), false));
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        for (leafOctant in leafOctants) {
            val octantAtoms = leafOctant.atoms;
            for (atoms in octantAtoms) {
                props.accumulate(atoms, this);
            }
            // print start positions
            printPositions(0, octantAtoms, startPosPrinter);
        }
        reduceAndPrintProperties(0, props);
    }

    private def printEndPositions(time:Double) {
        // print end positions
        val timeInt = (time) as Int;
        val endPosPrinter = new Printer(new FileWriter(new File("positions_" + timeInt + ".dat"), true));
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        for (leafOctant in leafOctants) {
            printPositions(time, leafOctant.atoms, endPosPrinter);
        }
    }

    private def printEnergies(time:Double) {
        // print end positions
        val timeInt = (time) as Int;
        val energiesPrinter = new Printer(new FileWriter(new File("energies_" + timeInt + ".dat"), true));
        var stats:Statistics = new Statistics();
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        for (leafOctant in leafOctants) {
            val s = printEnergies(time, leafOctant.atoms, energiesPrinter);
            stats.add(s);
        }
        energiesPrinter.printf("mean energy %.5g stddev %.4g", stats.mean(), stats.stdDev());
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param step the current step
     * @param dt time in s
     * @param current the current array at this place
     * @param whether to accumulate system properties for this step
     * @param props the system properties to evaluate
     */
    public def mdStepLocal(step:Int, dt:Double, current:Array[Double], accumProps:Boolean, props:SystemProperties, timer:StatisticalTimer) {
        timer.start(0);
        timer.start(2);
        fmm.reassignAtoms(step);
        timer.stop(2);

        timer.start(1);
        val pairwiseEnergy = fmm.calculateEnergyLocal();
        timer.stop(1);

        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        finish for (leafOctant in leafOctants) async {
            for (atom in leafOctant.atoms) {
                // timestep using Boris integrator
                val chargeMassRatio = atom.charge / atom.mass * CHARGE_MASS_FACTOR;

                // get the electric field at ion centre due to other ions
                /*
                var E:Vector3d = Vector3d.NULL;
                for (j in 0..(boxAtoms.size-1)) {
                    if (i == j) continue;
                    val atomJ = boxAtoms(j);
                    if (atomJ != null) {
                        val r = atom.centre - atomJ.centre;
                        val r2 = r.lengthSquared();
                        val absR = Math.sqrt(r2);
                        val atomContribution = COULOMB_FACTOR * atomJ.charge / r2 / absR * r;
                        //Console.OUT.println(j + " => " + i + " = " + atomContribution);
                        E = E + atomContribution;
                    }
                }
                var Ef:Vector3d = getElectrostaticField(atom.centre);
                //Console.OUT.println("Ef(" + i + ") = " + Ef.length() + " Ej = " + E.length() + " proportion " + Ef.length()/E.length());
                E += Ef;
                */
                val E = getElectrostaticField(atom.centre);

                //Console.OUT.println("E = " + E);
                val halfA = 0.5 * dt * chargeMassRatio * E;
                val vMinus = atom.velocity + halfA; // m

                //val larmorFreq = chargeMassRatio * magB;
                val t = 0.5 * dt * chargeMassRatio * B;
                val vPrime = vMinus + vMinus.cross(t);

                val magt2 = t.lengthSquared();
                val s = t * (2.0 / (1.0 + magt2));
                val vPlus = vMinus + vPrime.cross(s);

                atom.velocity = vPlus + halfA;
            }
        }

        Team.WORLD.barrier(here.id);

        var currentLocal:Double = 0.0;

        // TODO async
        //var placeAtoms:Int = 0;
        for (leafOctant in leafOctants) {
            val boxAtoms = leafOctant.atoms;
            //placeAtoms += boxAtoms.size;
            for (i in 0..(boxAtoms.size()-1)) {
                val atom = boxAtoms(i);
                atom.centre = atom.centre + atom.velocity * dt;
                //Console.OUT.print(atom.centre.i + " " + atom.centre.j + " " + atom.centre.k + " ");

                if (Math.abs(atom.centre.i) > edgeLength/2.0
                 || Math.abs(atom.centre.j) > edgeLength/2.0
                 || Math.abs(atom.centre.k) > edgeLength/2.0) {
                    // ion lost to wall
                    boxAtoms(i) = null;
                } else {
                    currentLocal += getImageCurrent(atom);
                    if (accumProps) {
                        props.accumulate(atom, this);
                    }
                }
            }
        }
        //Console.OUT.println("at " + here + " placeAtoms: " + placeAtoms);
        current(step) = currentLocal;
        timer.stop(0);
    }

    /** 
     * @return the position-dependent electrostatic field due to trapping plates
     * @see Guan & Marshall eq. 44
     */
    public @Inline def getElectrostaticField(p:Point3d):Vector3d {
        return Vector3d(p.i*eNorm, p.j*eNorm, -2.0*p.k*eNorm);
/*
        Han & Shin eq. 7-8
        var sumTerms:Double = 0.0;
        var sign:Double = -1.0;
        for (l in 0..10) {
            for (m in 0..10) {
                val k_lm = Math.sqrt((2*l+1)*(2*l+1) + (2*m+1)*(2*m+1));
                sumTerms += sign * Math.cosh(k_lm * Math.PI * z/a) / Math.cosh(k_lm * Math.PI / 2.0)
                     * Math.cos((2*l+1)*Math.PI*x/a)/(2*l+1)
                     * Math.cos((2*m+1)*Math.PI*y/a)/(2*m+1);
                sign *= -1.0;
            }
        }
        val V_0 = 0.0; // TODO side plate potential
        val phi = V_0 + 16.0*(V_T - V_0) / Math.PI*Math.PI * sumTerms;
        return phi;
*/
    }

    /** 
     * @return the position-dependent electrostatic potential due to trapping plates
     * @see Guan & Marshall eq. 19
     */
    public @Inline def getElectrostaticPotential(p:Point3d):Double {
        val potential = V_T * (/*1.0/3.0*/ -0.5 * eNorm * (p.i*p.i + p.j*p.j - 2.0*p.k*p.k));
        return potential;
    }

    /**
     * @return the image current induced by the given ion
     * @see Guan & Marshall eq. 69-70
     */
    public @Inline def getImageCurrent(ion:MMAtom):Double {
        return ion.charge * ion.velocity.j * eFieldNorm;
    }

    private def printPositions(time:Double, myAtoms:ArrayList[MMAtom], printer:Printer) {
        val timeInt = time as Int;
        for ([i] in 0..(myAtoms.size()-1)) {
            val atom = myAtoms(i);
            if (atom != null) {
                printer.printf("%i %i %s %12.8f %12.8f %12.8f\n", here.id, i, atom.symbol, atom.centre.i*1.0e3, atom.centre.j*1.0e3, atom.centre.k*1.0e3);
            }
        }
    }

    /** 
     * Prints the energies of all particles in the list.  Also returns
     * Statistics on the energy.
     */
    private def printEnergies(time:Double, myAtoms:ArrayList[MMAtom], printer:Printer):Statistics {
        val timeInt = time as Int;
        var stats:Statistics = new Statistics();
        for ([i] in 0..(myAtoms.size()-1)) {
            val atom = myAtoms(i);
            if (atom != null) {
                val energy = 0.5 * 1.66053892173e-12 /* Da->kg * 10^15fJ */ * atom.mass * atom.velocity.lengthSquared();
                printer.printf("%i %i %s %12.8f\n", here.id, i, atom.symbol, energy);
                stats.n++; 
                stats.sum += energy;
                stats.sumOfSquares += energy*energy;
            }
        }
        return stats;
    }

    private def reduceAndPrintProperties(time:Double, props:SystemProperties) {
        Team.WORLD.allreduce[Double](here.id, props.raw, 0, props.raw, 0, props.raw.size, Team.ADD);
        if (here == Place.FIRST_PLACE) {
            props.print(time);
        }
        props.reset();
    }


    private def printCurrent(timestep:Double, current:Rail[Double]) {
        val currentFilePrinter = new Printer(new FileWriter(new File("current.dat")));
        for ([i] in current) {
            current(i) *= 1.6021765314e-7; // e->C * 10^12; pA
            currentFilePrinter.printf("%16.8f\n", current(i));
            //currentFilePrinter.printf("%10.2f %16.8f\n", i*timestep, I);
        }
    }

    private def printMassSpectrum(timestep:Double, current:Rail[Double]) {
        val massSpectrum = new Array[Complex](current.size/2 + 1);
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(current.size, current, massSpectrum);
        FFTW.fftwExecute(plan);
        FFTW.fftwDestroyPlan(plan);

        val mzPrinter = new Printer(new FileWriter(new File("penning.mz")));

        val invN = 1.0e9/(current.size as Double);
        val sampleFreq = 1.0 / timestep * invN;

        val freq = new Array[Double](massSpectrum.size);
        val amplitude = new Array[Double](massSpectrum.size);
        mzPrinter.println("# Mass spectrum for FT-ICR.  Peaks at:");
        var peak:Double = 0.0;
        var increasing:Boolean = false;
        for ([i] in massSpectrum) {
            freq(i) = (i as Double) * sampleFreq;
            amplitude(i) = massSpectrum(i).abs();
            if (increasing && amplitude(i) <= peak) {
                increasing = false;
                if (i > 1) { // skip 'peak' at zero frequency
                    mzPrinter.printf("# %10.2f Hz\n", freq(i-1));
                }
            } else if (amplitude(i) > peak) {
                increasing = true;
            }
            peak = amplitude(i);
        }
        mzPrinter.println("#");
        mzPrinter.printf("#%9s %12s\n", "freq (Hz)", "amplitude");

        for (i in 0..(massSpectrum.size/4)) {
            mzPrinter.printf("%10.2f %.4g\n", freq(i), amplitude(i));
            //currentFilePrinter.printf("%10.2f %16.8f\n", i*timestep, I);
        }
    }

    static class SystemProperties { 
        public val raw:Rail[Double];
        public def this() {
            raw = new Array[Double](7);
        }

        /**
         * Sets all properties to zero.
         */
        public def reset() {
            raw.clear();
        }

        /**
         * Accumulates the properties for the given atom.
         */
        public @Inline def accumulate(atom:MMAtom, trap:PenningTrap) {
            raw(0)++;
            raw(1) += atom.centre.i;
            raw(2) += atom.centre.j;
            raw(3) += atom.centre.k;
            raw(4) += atom.mass * atom.velocity.lengthSquared();
            raw(5) += atom.charge * trap.getElectrostaticPotential(atom.centre);
            raw(6) += trap.getImageCurrent(atom);
        }

        /**
         * Prints the system properties.  Called at place 0 after reduction across all places.
         */
        public @Inline def print(time:Double) {
            val numAtoms = raw(0);
            val meanX = raw(1) / numAtoms * 1e3;
            val meanY = raw(2) / numAtoms * 1e3;
            val meanZ = raw(3) / numAtoms * 1e3;
            val Ek = raw(4) * 0.5 * 1.66053892173e-12; // Da->kg * 10^15 fJ
            val Ep = raw(5) * 1.6021765314e-4; // e->C * 10^15; fJ
            val E = Ek + Ep;
            val I = raw(6) * 1.6021765314e-7; // e->C * 10^12; pA

            Console.OUT.printf("%10.1f %8i %15.8f %15.8f %15.8f ", 
                time, 
                numAtoms as Int,
                meanX, 
                meanY, 
                meanZ);
            Console.OUT.printf("%15.8f %15.8f %15.8f %15.8f\n", 
                Ek, 
                Ep, 
                E, 
                I);
        }

        public def getCurrent() {
            return raw(6);
        }

        public static def printHeader() {
            Console.OUT.printf("%10s %8s %15s %15s %15s ",
                "time (ns)", 
                "num_ions",
                "mean_X (mm)",
                "mean_Y (mm)",
                "mean_Z (mm)");
            Console.OUT.printf("%15s %15s %15s %15s\n",
                "Ek (fJ)",
                "Ep (fJ)",
                "E (fJ)",
                "I (pA)"); 
        }
    }

    public static class Statistics {
        public var n:Int;
        public var sum:Double;
        public var sumOfSquares:Double;

        public def add(a:Statistics) {
            n += a.n;
            sum += a.sum;
            sumOfSquares += a.sumOfSquares;
        }

        public def mean() = sum / n;

        public def variance() {
            val mean = mean();
            return (sumOfSquares / n - mean*mean);
        }

        public def stdDev() {
            return Math.sqrt(variance());
        } 
    }
}

