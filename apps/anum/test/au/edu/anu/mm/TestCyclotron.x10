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
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Simulates an ion packet travelling in a circular path in a Penning trap
 * for the purposes of Fourier Transform Ion Cyclotron Resonance (FT-ICR) mass
 * spectroscopy.
 * @author milthorpe
 */
public class TestCyclotron extends TestElectrostatic {
    public static def main(args : Array[String](1)) {

        Console.OUT.println("Testing simple HF morse oscillator.");

        // equilibrium bond length is 0.09169nm; start with displacement of 0.01nm
        val hydrogen = new MMAtom("H", Point3d(0, 0.1, 0.0), 1.0);
        hydrogen.velocity = Vector3d(0.0, 10.0, 0.0);
        val atoms = new Array[MMAtom](1);
        atoms(0) = hydrogen;
        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val anumm = new Anumm(distAtoms, new PenningForceField());
        anumm.mdRun(0.2, 2000);
    }
}

