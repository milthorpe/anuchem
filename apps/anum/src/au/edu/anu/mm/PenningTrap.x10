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
 */
public class PenningTrap {
    public static val LENGTH_FACTOR = 1.0e-9;
    public static val CHARGE_MASS_FACTOR = 9.64853364e7; // conversion of q/m from e/Da to C/kg
    public static val MASS_CHARGE_FACTOR = 1.036426941e-8; // conversion of m/q from Da/e to kg/C
    static val ALPHA_PRIME = 2.77373; // geometric factor for a cubic trap (Guan and Marshall eq. 59)
    static val BETA_PRIME = 0.72167; // electric field constant for detection/excition (Guan and Marshall eq. 66)
    static val EDGE_LENGTH = 0.0508; // edge length l = 5.08
    static val COULOMB_FACTOR = 3.501233e-27; // 1 / (4 PI e0 * CHARGE_FACTOR^2)

    private var numAtoms:Int;

    /** The atoms in the simulation, divided up into a distributed array of Arrays, one for each place. */
    private val atoms:DistArray[Rail[MMAtom]](1);

    /** Normalization factor for electric field. */
    static val E_NORM = 1255.65; // ALPHA_PRIME / (EDGE_LENGTH^2)

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
     * @param properties system properties to print at each timestep
     */
    public def this(numAtoms:Int,
                    atoms:DistArray[Rail[MMAtom]](1),
                    trappingPotential:Double,
                    magneticField:Vector3d) {
        this.numAtoms = numAtoms;
        this.atoms = atoms;
        this.V_T = trappingPotential;
        this.B = magneticField;
        this.magB = magneticField.magnitude();
    }

    public def getAtoms() = atoms;

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in ns
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Long, logSteps:Long) {
        Console.OUT.println("# Timestep = " + timestep + "ns, number of steps = " + numSteps);
        val timer = new Timer(2);
 
        SystemProperties.printHeader();
           
        finish ateach(placeId in atoms) {
            var step : Long = 0L;
            val myAtoms = atoms(placeId);
            val props = new SystemProperties();
            for (i in 0..(myAtoms.size-1)) {
                if (myAtoms(i) != null) {
                    props.accumulate(myAtoms(i), this);
                }
            }
            if (here == Place.FIRST_PLACE) {
                props.print(0, numAtoms);
            }
            while(step < numSteps) {
                step++;
                props.reset();
                mdStep(timestep, myAtoms, props);
                if (step % logSteps == 0L) {
                    Team.WORLD.allreduce[Double](here.id, props.raw, 0, props.raw, 0, props.raw.size, Team.ADD);
                    if (here == Place.FIRST_PLACE) {
                        props.print(timestep * step, numAtoms);
                    }
                }
            }
            printPositions(timestep * step, myAtoms);
        }
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param dt time in ns
     * @param myAtoms the atoms that are stored here
     * @param props the system properties to evaluate
     */
    public def mdStep(dt:Double, myAtoms:Rail[MMAtom], props:SystemProperties) {
        for (i in 0..(myAtoms.size-1)) {
            val atom = myAtoms(i);
            if (atom != null) {
    
                // timestep using Boris integrator
                val chargeMassRatio = atom.charge / atom.mass * CHARGE_MASS_FACTOR;

                var E:Vector3d = getElectrostaticField(atom.centre);
                /*
                //Console.OUT.println("E(" + i + ") = " + E);
                for (j in 0..(myAtoms.size-1)) {
                    if (i == j) continue;
                    val atomJ = myAtoms(j);
                    if (atomJ != null) {
                        val r = atom.centre - atomJ.centre;
                        val r2 = r.lengthSquared();
                        val absR = Math.sqrt(r2);
                        val atomContribution = COULOMB_FACTOR * atom.charge * atomJ.charge / r2 / absR * r;
                        //Console.OUT.println(j + " => " + i + " = " + atomContribution);
                        E = E + atomContribution;
                    }
                }
                */

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

        // TODO global barrier
        for (i in 0..(myAtoms.size-1)) {
            val atom = myAtoms(i);
            if (atom != null) {
                atom.centre = atom.centre + atom.velocity * dt * 1.0e-9;
                //Console.OUT.print(atom.centre.i + " " + atom.centre.j + " " + atom.centre.k + " ");

                if (Math.abs(atom.centre.i) > EDGE_LENGTH/2.0
                 || Math.abs(atom.centre.j) > EDGE_LENGTH/2.0
                 || Math.abs(atom.centre.k) > EDGE_LENGTH/2.0) {
                    // ion lost to wall
                    myAtoms(i) = null;
                    numAtoms--;
                }

                props.accumulate(atom, this);
            }
        }
    }

    /** 
     * @return the position-dependent electrostatic field due to trapping plates
     * @see Guan & Marshall eq. 44
     */
    public static @Inline def getElectrostaticField(p:Point3d):Vector3d {
        return Vector3d(p.i*E_NORM, p.j*E_NORM, -2.0*p.k*E_NORM);
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
        val potential = V_T * (/*1.0/3.0*/ -0.5 * E_NORM * (p.i*p.i + p.j*p.j - 2.0*p.k*p.k));
        return potential;
    }

    /**
     * @return the image current induced by the given ion
     * @see Guan & Marshall eq. 69-70
     */
    public @Inline static def getImageCurrent(ion:MMAtom):Double {
        val eImage = -BETA_PRIME / EDGE_LENGTH;
        return ion.charge * ion.velocity.j * eImage;
    }

    private def printPositions(time:Double, myAtoms:Rail[MMAtom]) {
        val timeInt = time as Int;
        val posFilePrinter = new Printer(new FileWriter(new File("positions_" + timeInt + ".dat")));
        for ([i] in myAtoms) {
            val atom = myAtoms(i);
            if (atom != null) {
                posFilePrinter.printf("%i %i %s %12.8f %12.8f %12.8f\n", here.id, i, atom.symbol, atom.centre.i*1.0e3, atom.centre.j*1.0e3, atom.centre.k*1.0e3);
            }
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
            raw(5) += getImageCurrent(atom);
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

            Console.OUT.printf("%10.2f %8i", 
                time, 
                numAtoms);
            Console.OUT.printf("%16.8f %16.8f %16.8f ", 
                meanX, 
                meanY, 
                meanZ);
            Console.OUT.printf("%16.8f %16.8f %16.8f %16.8f\n", 
                Ek, 
                Ep, 
                E, 
                I);
        }

        public static def printHeader() {
            Console.OUT.printf("%10s %8s", 
                "ns", 
                "num_ions");
            Console.OUT.printf("%16s %16s %16s ",
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

