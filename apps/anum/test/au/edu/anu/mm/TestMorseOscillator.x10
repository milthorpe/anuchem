/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */
package au.edu.anu.mm;

import x10.regionarray.Dist;
import x10.regionarray.DistArray;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;

/**
 * Tests ANU Molecular Mechanics with a simple Morse oscillator: HF.
 * A test problem borrowed from Chapter 1 of
 * Herman J.C. Berendsen, "Simulating the Physical World"
 * Cambridge University Press, 2007
 */
public class TestMorseOscillator {
    public static def main(args:Rail[String]) {
        Console.OUT.println("Testing simple HF morse oscillator.");

        // equilibrium bond length is 0.09169nm; start with displacement of 0.01nm
        val hydrogen = new MMAtom(DiatomicForceField.SPECIES_H, Point3d(0.10169, 0.0, 0.0), 1.0079, 0.0);
        val fluorine = new MMAtom(DiatomicForceField.SPECIES_F, Point3d(0.0, 0.0, 0.0), 18.9984, 0.0);
        val atoms = new Rail[MMAtom](2);
        atoms(0) = hydrogen;
        atoms(1) = fluorine;
        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeUnique());
        distAtoms(0) = atoms;

        val diatomicPotentials = new Rail[DiatomicPotential](1L, 
           (Long) => new DiatomicMorsePotential(hydrogen, fluorine, 0.09169, 569.87));
        val anumm = new Anumm(distAtoms, new DiatomicForceField(diatomicPotentials), true);
        anumm.mdRun(0.2, 200);
    }
}

