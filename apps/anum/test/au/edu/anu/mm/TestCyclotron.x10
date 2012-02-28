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

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.mm.SystemProperties;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;

/**
 * Simulates an ion packet travelling in a circular path in a Penning trap
 * for the purposes of Fourier Transform Ion Cyclotron Resonance (FT-ICR) mass
 * spectroscopy.
 * @author milthorpe
 */
public class TestCyclotron extends TestElectrostatic {
    public static def main(args : Array[String](1)) {
        var v:Double = 1.0;
        var timesteps:Int = 200000;
        if (args.size > 0) {
            v = Double.parseDouble(args(0));
            if (args.size > 1) {
                timesteps = Int.parseInt(args(1));
            }
        }

        Console.OUT.println("Testing cyclotron: initial velocity: + " + v + "");

        // start with displacement of 0.01nm
        val hydrogen = new MMAtom("H", Point3d(0, 0.00, 0.001), 1.0);
        hydrogen.velocity = Vector3d(0.0, v, 0.0);
        val atoms = new Array[MMAtom](1);
        atoms(0) = hydrogen;
        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val onePFs = [OneParticleFunction("mean_Z", (a:MMAtom) => a.centre.k)];
        val anumm = new Anumm(distAtoms, new PenningForceField(new Vector3d(-1.0, 0.0, 0.0)), new SystemProperties(onePFs));
        anumm.mdRun(0.2, timesteps, 100);
    }
}

