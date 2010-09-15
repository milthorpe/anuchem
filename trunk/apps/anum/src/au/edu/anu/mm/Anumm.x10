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
package au.edu.anu.mm;

import x10.util.*;
import x10.io.IOException;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.Molecule;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;
import au.edu.anu.mm.uff.UniversalForceField;

/**
 * ANU-Mechanics, a molecular mechanics simulator.
 * Performs molecular mechanics using the velocity-Verlet algorithm.
 */
public class Anumm {
    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : DistArray[ValRail[MMAtom]](1);

    /** The force field applied to the atoms in this simulation. */
    private global val forceField : ForceField;

    public def this(atoms: DistArray[ValRail[MMAtom]](1),
                    forceField: ForceField) {
        this.atoms = atoms;
        this.forceField = forceField;
    }

    public def getAtoms() = atoms;

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in fs (=ps/1000)
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep : Double, numSteps : Long) {
        Console.OUT.println("# Timestep = " + timestep + "fs, number of steps = " + numSteps);
        Console.OUT.println("0.0 ");
        forceField.getPotentialAndForces(atoms); // get initial forces
        var steps : Long = 0;
        while(steps < numSteps) {
            steps++;
            Console.OUT.print("\n" + timestep * steps + " ");
            mdStep(timestep);
        }
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param timestep time in fs (=ps/1000)
     */
    public def mdStep(timestep : Double) {
        val t = timestep * 0.001;
        finish ateach([p] in atoms) {
            val myAtoms = atoms(p);
            for([i] in 0..myAtoms.length-1) {
                val atom = myAtoms(i);
                val invMass = 1.0 / forceField.getAtomMass(atom.symbol);
                atom.velocity = atom.velocity + 0.5 * t * invMass * atom.force;
                atom.centre = atom.centre + atom.velocity * t;
            }
        }
        forceField.getPotentialAndForces(atoms);
        finish ateach([p] in atoms) {
            val myAtoms = atoms(p);
            for([i] in 0..myAtoms.length-1) {
                val atom = myAtoms(i);
                val invMass = 1.0 / forceField.getAtomMass(atom.symbol);
                atom.velocity = atom.velocity + 0.5 * t * invMass * atom.force;
                Console.OUT.print(atom + " ");
            }
        }
    }

    public static def main(args : Array[String](1)) {
        var structureFileName : String = null;
        var timestep : Double = 0.2;
        var numSteps : Int = 200;
        if (args.length > 0) {
            structureFileName = args(0);
            if (args.length > 1) {
                timestep = Double.parseDouble(args(1));
                if (args.length > 2) {
                    numSteps = Int.parseInt(args(2));
                }
            }
        } else {
            Console.ERR.println("usage: anumm structureFile [timestep(fs)] [numSteps]");
            return;
        }
        var moleculeTemp : Molecule[MMAtom] = null;
        if (structureFileName.length() > 4) {
            val fileExtension = structureFileName.substring(structureFileName.length()-4, structureFileName.length());
            if (fileExtension.equals(".xyz")) {
                try {
                    moleculeTemp = new XYZStructureFileReader(structureFileName).readMolecule();
                } catch (e:IOException) {
                    Console.ERR.println(e);
                }
            } else if (fileExtension.equals(".mol") || fileExtension.equals(".sdf")) {
                try {
                    moleculeTemp = new MOLStructureFileReader(structureFileName).readMolecule();
                } catch (e:IOException) {
                    Console.ERR.println(e);
                }
            }
        }
        if (moleculeTemp == null) {
            Console.ERR.println("error: could not read structure file: " + structureFileName);
            return;
        }
        val molecule = moleculeTemp;
        Console.OUT.println("# MD for " + molecule.getName() + ": " + molecule.getAtoms().size() + " atoms");
        val anumm = new Anumm(assignAtoms(molecule), new UniversalForceField());
        anumm.mdRun(timestep, numSteps);

        return;
    }

    /**
     * Partitions the molecule amongst all places, returning a distributed
     * array of ValRail[MMAtom], one ValRail for each place.  
     * MD requires that the atoms have already been distributed. 
     */
    public static def assignAtoms(molecule : Molecule[MMAtom]) : DistArray[ValRail[MMAtom]](1) {
        val tempAtoms = DistArray.make[GrowableRail[MMAtom]](Dist.makeUnique(Place.places), (Point) => new GrowableRail[MMAtom]());
        val atomList = molecule.getAtoms();
        val maxExtent = molecule.getMaxExtent();
        finish for (var i : Int = 0; i < atomList.size(); i++) {
            val atom = atomList(i);
            val p = getPlaceId(atom.centre.i, atom.centre.j, atom.centre.k, maxExtent);
            Console.OUT.println(atom + " to " + p);
            async (Place.places(p)) {
                val remoteAtom = new MMAtom(atom);
                atomic { (tempAtoms(p) as GrowableRail[MMAtom]!).add(remoteAtom); }
            }
        }
        val atoms = DistArray.make[ValRail[MMAtom]](Dist.makeUnique(Place.places), ((p) : Point) => (tempAtoms(p) as GrowableRail[MMAtom]!).toValRail());
        return atoms;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    private static safe def getPlaceId(x : Double, y : Double, z : Double, size : Double) : Int {
        return ((x / (size * 2) + 0.5) * Place.MAX_PLACES) as Int;
    }

}

