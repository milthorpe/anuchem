package au.edu.anu.mm;

import x10.util.*;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * ANU-Mechanics, a molecular mechanics simulator.
 * Performs molecular mechanics using the velocity-Verlet algorithm.
 */
public class Anum {
    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : DistArray[ValRail[MMAtom]](1);

    /** The force field applied to the atoms in this simulation. */
    private global val forceField : ForceField;

    public def this(atoms: DistArray[ValRail[MMAtom]](1),
                    forceField: ForceField) {
        this.atoms = atoms;
        this.forceField = forceField;
    }

    /**
     * Perform a molecular mechanics run on the system of atoms
     * for the given number and length of timesteps.
     * @param timestep length in fs (=ps/1000)
     * @param numSteps number of timesteps to simulate
     */
    public def mdRun(timestep : Double, numSteps : Long) {
        forceField.getPotentialAndForces(atoms); // get initial forces
        var steps : Long = 0;
        while(steps < numSteps) {
            Console.OUT.print("\n" + timestep * steps + " ");
            mdStep(timestep);
            steps++;
        }
    }

    /**
     * Performs a single molecular dynamics timestep
     * using the velocity-Verlet algorithm. 
     * @param timestep time in fs (=ps/1000)
     */
    public def mdStep(timestep : Double) {
        val t = timestep * 0.001;
        finish ateach((p) in atoms) {
            val myAtoms = atoms(p);
            for((i) in 0..myAtoms.length-1) {
                val atom = myAtoms(i);
                val invMass = 1.0 / atom.mass;
                atom.velocity = atom.velocity + 0.5 * t * invMass * atom.force;
                atom.centre = atom.centre + atom.velocity * t;
            }
        }
        forceField.getPotentialAndForces(atoms);
        finish ateach((p) in atoms) {
            val myAtoms = atoms(p);
            for((i) in 0..myAtoms.length-1) {
                val atom = myAtoms(i);
                val invMass = 1.0 / atom.mass;
                atom.velocity = atom.velocity + 0.5 * t * invMass * atom.force;
                //Console.OUT.print(atom + " ");
            }
        }
    }
}

