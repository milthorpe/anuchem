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

import x10.util.ArrayList;
import x10.util.Team;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * A simulation of charged ions in cyclotron motion in a Penning trap.
 * Equations of motion are integrated using the Boris scheme.
 */
public class PenningTrap {
    static val CHARGE_MASS_FACTOR = 9.64853105e-8; // conversion of q/m from e/Da to C/kg / 1e15

    private val numAtoms:Int;

    /** The atoms in the simulation, divided up into a distributed array of Arrays, one for each place. */
    private val atoms:DistArray[Rail[MMAtom]](1);

    /** The static homogeneous magnetic field B, in Teslas. */
    public val B:Vector3d;
    /** The scalar magnitude of B. */
    private val magB:Double;

    // static E = ??

    /** The system properties to be calculated at each log timestep. */
    private val properties:SystemProperties;

    public def this(numAtoms:Int,
                    atoms:DistArray[Rail[MMAtom]](1),
                    magneticField:Vector3d,
                    properties:SystemProperties) {
        this.numAtoms = numAtoms;
        this.atoms = atoms;
        this.B = magneticField;
        this.magB = magneticField.magnitude();
        this.properties = properties;
    }

    public def getAtoms() = atoms;

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in fs
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep:Double, numSteps:Long, logSteps:Long) {
        Console.OUT.println("# Timestep = " + timestep + "fs, number of steps = " + numSteps);

        Console.OUT.printf("%12s ", "us");
        val funcs = properties.oneParticleFunctions;
        for (i in 0..(funcs.size-1)) {
            Console.OUT.printf("%16s ", funcs(i).first);
        }
        Console.OUT.println();

        finish ateach(placeId in atoms) {
            var step : Long = 0;
            val myAtoms = atoms(placeId);
            printProperties(timestep, step, myAtoms);
            while(step < numSteps) {
                step++;
                mdStep(timestep, myAtoms);
                if (step % logSteps == 0L) {
                    printProperties(timestep, step, myAtoms);
                }
            }
        }
    }

    /**
     * Print current system properties.
     * @param timestep length in fs
     * @param currentStep current time in number of steps from start
     * @param myAtoms for which to calculate properties
     */
    private def printProperties(timestep:Double, currentStep:Long, myAtoms:Rail[MMAtom]) {
        val propertySums = properties.calculateExpectationValues(myAtoms);
        
        Team.WORLD.allreduce[Double](here.id, propertySums, 0, propertySums, 0, propertySums.size, Team.ADD);
        for (i in 0..(propertySums.size-1)) {
            propertySums(i) /= (numAtoms as Double);
        }
        if (here == Place.FIRST_PLACE) {
            Console.OUT.printf("%12.6f ", timestep * currentStep / 1e9);
            for (i in 0..(propertySums.size-1)) {
                Console.OUT.printf("%16.8f ", propertySums(i));
            }
            Console.OUT.println();
        }
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param dt time in fs
     */
    public def mdStep(dt:Double, myAtoms:Rail[MMAtom]) {
        for (i in 0..(myAtoms.size-1)) {
            val atom = myAtoms(i);
    
            // timestep using Boris integrator
            val invMass = 1.0 / getAtomMass(atom.symbol);

            //val halfA = 0.5 * dt * atom.charge * E * invMass;
            val vMinus = atom.velocity; // + halfA;

            val chargeMassRatio = atom.charge * invMass * CHARGE_MASS_FACTOR;

            //val larmorFreq = chargeMassRatio * magB;
            val t = dt * 0.5 * chargeMassRatio * B;
            val vPrime = vMinus + vMinus.cross(t);

            val magt2 = t.lengthSquared();
            val s = t * (2.0 / (1.0 + magt2));
            val vPlus = vMinus + vPrime.cross(s);

            atom.velocity = vPlus; // + halfA;

            atom.centre = atom.centre + atom.velocity*dt;
            //Console.OUT.print(atom.centre.i + " " + atom.centre.j + " " + atom.centre.k + " ");
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

    /** @return atomic/molar mass in atomic units */
    public static def getAtomMass(symbol : String) : Double {
        if (symbol.equals("H")) {
            return 1.0079;
        } else if (symbol.equals("F")) {
            return 18.9984;
        } else if (symbol.equals("CH3CO")) {
            // acetaldehyde cation
            return 43.04462;
        } else if (symbol.equals("HCO")) {
            // HCO+ cation 
            return 29.0182;
        } else {
            throw new IllegalArgumentException("no atom mass found for symbol " + symbol);
        }
    }

}

