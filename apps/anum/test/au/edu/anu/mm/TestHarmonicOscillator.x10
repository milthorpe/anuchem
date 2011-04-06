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

import x10.util.Pair;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * Tests ANU Molecular Mechanics with a simple Morse oscillator: HF.
 * A test problem borrowed from Chapter 1 of
 * Herman J.C. Berendsen, "Simulating the Physical World"
 * Cambridge University Press, 2007
 * @author milthorpe
 */
public class TestHarmonicOscillator {
    public static def main(args : Array[String](1)) {

        Console.OUT.println("Testing simple HF harmonic oscillator.");

        // equilibrium bond length is 0.09169nm; start with displacement of 0.01nm
        val hydrogen = new MMAtom("H", Point3d(0.10169, 0.0, 0.0), 0.0);
        val fluorine = new MMAtom("F", Point3d(0.0, 0.0, 0.0), 0.0);
        val atoms = new Array[MMAtom](2);
        atoms(0) = hydrogen;
        atoms(1) = fluorine;
        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val diatomicPotentials = new Array[DiatomicPotential](1, 
           (Int) => new DiatomicHarmonicPotential(hydrogen, fluorine, 0.09169, 582000));
        val anumm = new Anumm(distAtoms, new DiatomicForceField(diatomicPotentials));
        anumm.mdRun(0.2, 200);
    }
}

