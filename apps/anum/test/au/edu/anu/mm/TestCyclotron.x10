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

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.mm.SystemProperties;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;

/**
 * Simulates an ion packet travelling in a circular path in a Penning trap
 * for the purposes of Fourier Transform Ion Cyclotron Resonance (FT-ICR) mass
 * spectroscopy.
 * Parameters to match simulation reported in:
 * @see Seung-Jin Han and Seung Koo Shin (1997)
 *  "Space-Charge Effects and Fourier transform ion cyclotron resonance signals:
 *   experimental observations and three-dimensional trajectory simulations"
 *   J. Am. Soc. Mass Spectrometry 8 (4), 319-326 
 * @author milthorpe
 */
public class TestCyclotron extends TestElectrostatic {
    public static def main(args : Array[String](1)) {
        var B:Double = 0.7646;
        var dt:Double = 50000.0; // timestep in fs
        var V:Double = 1.0;
        var timesteps:Int = 80000; // number of timesteps
        var logSteps:Int = 100;
        if (args.size > 0) {
            B = Double.parseDouble(args(0));
            if (args.size > 1) {
                dt = Double.parseDouble(args(1));
                if (args.size > 2) {
                    V = Double.parseDouble(args(2));
                    if (args.size > 3) {
                        timesteps = Int.parseInt(args(3));
                        if (args.size > 4) {
                            logSteps = Int.parseInt(args(4));
                        }
                    }
                }
            }
        }

        Console.OUT.println("# Testing cyclotron: trapping potential: + " + V + " magnetic field: " + B);

        val r = new Random();
        val acetaldehyde = new MMAtom("CH3CO", Point3d(0.0, 0.0, 0.0), 1.0);
        acetaldehyde.velocity = Vector3d(r.nextDouble()*10, r.nextDouble()*10, 0.0);
        val atoms = new Array[MMAtom](1);
        atoms(0) = acetaldehyde;
        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val onePFs = [OneParticleFunction("mean_X", (a:MMAtom) => a.centre.i), OneParticleFunction("mean_Y", (a:MMAtom) => a.centre.j), OneParticleFunction("Ek", (a:MMAtom) => PenningTrap.getAtomMass(a.symbol) * a.velocity.lengthSquared())];
        val trap = new PenningTrap(1, distAtoms, V, new Vector3d(0.0, 0.0, B), new SystemProperties(1, onePFs));
        trap.mdRun(dt, timesteps, logSteps);
    }
}

