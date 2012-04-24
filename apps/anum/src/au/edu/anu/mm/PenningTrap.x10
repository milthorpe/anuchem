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
import au.edu.anu.util.Timer;
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

    val fmm:Fmm3d;

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
                    fmmDensity:Double,
                    fmmNumTerms:Int) {
        this.numAtoms = numAtoms;
        val trapRadius = edgeLength / 2.0;
        this.fmm = new Fmm3d(fmmDensity, fmmNumTerms, 1, Point3d(-trapRadius, -trapRadius, -trapRadius), edgeLength, numAtoms);
        fmm.assignAtomsToBoxes(atoms);
        this.V_T = trappingPotential;
        this.B = magneticField;
        this.edgeLength = edgeLength;
        this.eNorm = V_T * PenningTrap.ALPHA_PRIME / (edgeLength * edgeLength);
        this.magB = magneticField.magnitude();
    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in ns
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Int, logSteps:Int) {
        Console.OUT.println("# Timestep = " + timestep + "ns, number of steps = " + numSteps + " logging every " + logSteps);
        val timer = new Timer(2);

        val current = new Array[Double](numSteps);
 
        SystemProperties.printHeader();

        val fmmBoxes = fmm.boxes(fmm.numLevels);
        val posPrinter = new Printer(new FileWriter(new File("positions_0.dat"), false));
        posPrinter.println("# positions at time 0");
        finish ateach(place in Dist.makeUnique()) {
            var step:Int = 0;
            val props = new SystemProperties();
            for ([x,y,z] in fmmBoxes.dist(here)) {
                val box = fmmBoxes(x,y,z) as FmmLeafBox;
                if (box != null) {
                    val myAtoms = box.getAtoms();
                    for (i in 0..(myAtoms.size-1)) {
                        if (myAtoms(i) != null) {
                            props.accumulate(myAtoms(i), this);
                        }
                    }
                    // print start positions
                    printPositions(timestep * step, box.getAtoms());
                }
            }
            Team.WORLD.allreduce[Double](here.id, props.raw, 0, props.raw, 0, props.raw.size, Team.ADD);
            if (here == Place.FIRST_PLACE) {
                props.print(0, numAtoms);
            }

            while(step < numSteps) {
                step++;
                props.reset();
                mdStepLocal(timestep, fmmBoxes, props);
                if (step % logSteps == 0) {
                    Team.WORLD.allreduce[Double](here.id, props.raw, 0, props.raw, 0, props.raw.size, Team.ADD);
                    if (here == Place.FIRST_PLACE) {
                        props.print(timestep * step, numAtoms);
                    }
                }
                if (here == Place.FIRST_PLACE) {
                    current(step) = props.getCurrent();
                }
            }
            if (here == Place.FIRST_PLACE) {
                printCurrent(timestep, current);
                printMassSpectrum(timestep, current);
            }

            // print end positions
            val timeInt = (timestep * step) as Int;
            val endPosPrinter = new Printer(new FileWriter(new File("positions_" + timeInt + ".dat"), false));
            endPosPrinter.println("# positions at time " + timeInt);
            for ([x,y,z] in fmmBoxes.dist(here)) {
                val box = fmmBoxes(x,y,z) as FmmLeafBox;
                if (box != null) {
                    printPositions(timestep * step, box.getAtoms());
                }
            }
        }

    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param dt time in ns
     * @param fmmBoxes the lowest level boxes in the FMM tree
     * @param props the system properties to evaluate
     */
    public def mdStepLocal(dt:Double, fmmBoxes:DistArray[FmmBox](3), props:SystemProperties) {
        finish for ([x,y,z] in fmmBoxes.dist(here)) async {
            val box = fmmBoxes(x,y,z) as FmmLeafBox;
            if (box != null) {
                val myAtoms = box.getAtoms();
                for (i in 0..(myAtoms.size-1)) {
                    val atom = myAtoms(i);
                    if (atom != null) {
            
                        // timestep using Boris integrator
                        val chargeMassRatio = atom.charge / atom.mass * CHARGE_MASS_FACTOR;

                        // get the electric field at ion centre due to other ions
                        /*
                        var E:Vector3d = Vector3d.NULL;
                        for (j in 0..(myAtoms.size-1)) {
                            if (i == j) continue;
                            val atomJ = myAtoms(j);
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
                        val halfA = 0.5 * dt * 1.0e-9 * chargeMassRatio * E;
                        val vMinus = atom.velocity + halfA; // m

                        //val larmorFreq = chargeMassRatio * magB;
                        val t = 0.5 * dt * 1.0e-9 * chargeMassRatio * B;
                        val vPrime = vMinus + vMinus.cross(t);

                        val magt2 = t.lengthSquared();
                        val s = t * (2.0 / (1.0 + magt2));
                        val vPlus = vMinus + vPrime.cross(s);

                        atom.velocity = vPlus + halfA;
                    }
                }
            }
        }

        Team.WORLD.barrier(here.id);

        // TODO async
        for ([x,y,z] in fmmBoxes.dist(here)) {
            val box = fmmBoxes(x,y,z) as FmmLeafBox;
            if (box != null) {
                val myAtoms = box.getAtoms();
                var notNull:Int = myAtoms.size;
                for (i in 0..(myAtoms.size-1)) {
                    val atom = myAtoms(i);
                    if (atom != null) {
                        atom.centre = atom.centre + atom.velocity * dt * 1.0e-9;
                        //Console.OUT.print(atom.centre.i + " " + atom.centre.j + " " + atom.centre.k + " ");

                        if (Math.abs(atom.centre.i) > edgeLength/2.0
                         || Math.abs(atom.centre.j) > edgeLength/2.0
                         || Math.abs(atom.centre.k) > edgeLength/2.0) {
                            // ion lost to wall
                            myAtoms(i) = null;
                            numAtoms--;
                            notNull--;
                        } else {
                            val boxIndex = Fmm3d.getLowestLevelBoxIndex(atom.centre, fmm.lowestLevelDim, fmm.size);
                            if (boxIndex(0) != x || boxIndex(1) != y || boxIndex(2) != z) {
                                //Console.OUT.println("moving atom " + atom.centre + " from " + x+","+y+","+z + " to " + boxIndex);
                                at(fmmBoxes.dist(boxIndex)) {
                                    var destBox:FmmLeafBox = fmmBoxes(boxIndex) as FmmLeafBox;
                                    if (destBox == null) {
                                        destBox = new FmmLeafBox(fmm.numLevels, boxIndex(0), boxIndex(1), boxIndex(2), fmm.numTerms, Fmm3d.getParentForChild(fmm.boxes, fmm.numLevels, fmm.topLevel, x,y,z));
                                        val newAtoms = new Rail[MMAtom](1);
                                        newAtoms(0) = atom;
                                        destBox.setAtoms(newAtoms);
                                        fmmBoxes(boxIndex) = destBox;
                                    } else {
                                        val oldAtoms = destBox.getAtoms();
                                        val newAtoms = new Rail[MMAtom](oldAtoms.size+1);
                                        for (j in oldAtoms) {
                                            newAtoms(j) = oldAtoms(j);
                                        }
                                        newAtoms(newAtoms.size-1) = atom;
                                        destBox.setAtoms(newAtoms);
                                    }
                                }
                                myAtoms(i) = null;
                                notNull--;
                            }
                        }

                        props.accumulate(atom, this);
                    } else {
                        notNull--;
                    }
                }
                if (notNull == 0) {
                    //Console.OUT.println("deleting box " + x+","+y+","+z);
                    fmmBoxes(x,y,z) = null;
                } else if (notNull < myAtoms.size) {
                    // resize this box
                    //Console.OUT.println("resizing box " + x+","+y+","+z + " from " + myAtoms.size + " to " + notNull);
                    val newAtoms = new Rail[MMAtom](notNull);
                    var j:Int=0;
                    for (i in myAtoms) {
                        if (myAtoms(i) != null) {
                            newAtoms(j++) = myAtoms(i);
                        }
                    }
                    box.setAtoms(newAtoms);
                }
            }
        }
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
        val eImage = -BETA_PRIME / edgeLength;
        return ion.charge * ion.velocity.j * eImage;
    }

    private def printPositions(time:Double, myAtoms:Rail[MMAtom]) {
        val timeInt = time as Int;
        val posFilePrinter = new Printer(new FileWriter(new File("positions_" + timeInt + ".dat"), true));
        for ([i] in myAtoms) {
            val atom = myAtoms(i);
            if (atom != null) {
                posFilePrinter.printf("%i %i %s %12.8f %12.8f %12.8f\n", here.id, i, atom.symbol, atom.centre.i*1.0e3, atom.centre.j*1.0e3, atom.centre.k*1.0e3);
            }
        }
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

    /**
     * Partitions the molecule amongst all places, returning a distributed
     * array of Array[MMAtom], one Array for each place.  
     * MD requires that the atoms have already been distributed. 
     */
    public static def assignAtoms(molecule : Molecule[MMAtom]) : DistArray[Rail[MMAtom]](1) {
        val tempAtoms = DistArray.make[ArrayList[MMAtom]](Dist.makeUnique(), (Point) => new ArrayList[MMAtom]());
        val atomList = molecule.getAtoms();
        val maxExtent = molecule.getMaxExtent();
        finish for (var i : Int = 0; i < atomList.size(); i++) {
            val atom = atomList(i);
            val p = getPlaceId(atom.centre.i, atom.centre.j, atom.centre.k, maxExtent);
            //Console.OUT.println(atom + " to " + p);
            at(Place.place(p)) async {
                val remoteAtom = new MMAtom(atom);
                atomic { tempAtoms(p).add(remoteAtom); }
            }
        }
        val atoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique(), ([p] : Point) => tempAtoms(p).toArray());
        return atoms;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    private static def getPlaceId(x : Double, y : Double, z : Double, size : Double) : Int {
        return ((x / (size * 2) + 0.5) * Place.MAX_PLACES) as Int;
    }

    static class SystemProperties { 
        public val raw:Rail[Double];
        public def this() {
            raw = new Array[Double](6);
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
            raw(0) += atom.centre.i;
            raw(1) += atom.centre.j;
            raw(2) += atom.centre.k;
            raw(3) += atom.mass * atom.velocity.lengthSquared();
            raw(4) += atom.charge * trap.getElectrostaticPotential(atom.centre);
            raw(5) += trap.getImageCurrent(atom);
        }

        /**
         * Prints the system properties.  Called at place 0 after reduction across all places.
         */
        public @Inline def print(time:Double, numAtoms:Int) {
            val meanX = raw(0) / numAtoms * 1e3;
            val meanY = raw(1) / numAtoms * 1e3;
            val meanZ = raw(2) / numAtoms * 1e3;
            val Ek = raw(3) * 0.5 * 1.66053892173e-12; // Da->kg * 10^15 fJ
            val Ep = raw(4) * 1.6021765314e-4; // e->C * 10^15; fJ
            val E = Ek + Ep;
            val I = raw(5) * 1.6021765314e-7; // e->C * 10^12; pA

            Console.OUT.printf("%10.2f %8i %16.8f %16.8f %16.8f ", 
                time, 
                numAtoms,
                meanX, 
                meanY, 
                meanZ);
            Console.OUT.printf("%16.8f %16.8f %16.8f %16.8f\n", 
                Ek, 
                Ep, 
                E, 
                I);
        }

        public def getCurrent() {
            return raw(5);
        }

        public static def printHeader() {
            Console.OUT.printf("%10s %8s %16s %16s %16s ",
                "time (ns)", 
                "num_ions",
                "mean_X (mm)",
                "mean_Y (mm)",
                "mean_Z (mm)");
            Console.OUT.printf("%16s %16s %16s %16s\n",
                "Ek (fJ)",
                "Ep (fJ)",
                "E (fJ)",
                "I (pA)"); 
        }
    }
}

