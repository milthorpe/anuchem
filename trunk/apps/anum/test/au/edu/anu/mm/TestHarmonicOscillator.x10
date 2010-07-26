/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
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
    public static def main(args : Rail[String]!) {

        Console.OUT.println("Testing simple HF harmonic oscillator.");

        // equilibrium bond length is 0.09169nm; start with displacement of 0.01nm
        val hydrogen = new MMAtom("H", Point3d(0.10169, 0.0, 0.0), 1.0079, 0.0);
        val fluorine = new MMAtom("F", Point3d(0.0, 0.0, 0.0), 18.9984, 0.0);
        val atoms = Rail.make[MMAtom](2);
        atoms(0) = hydrogen;
        atoms(1) = fluorine;
        val distAtoms = DistArray.make[ValRail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms as ValRail[MMAtom];

        val diatomicPotentials : ValRail[DiatomicPotential] = ValRail.make[DiatomicPotential](1, 
           (Int) => new DiatomicHarmonicPotential(hydrogen, fluorine, 0.09169, 582000));
        val anumm = new Anumm(distAtoms, new DiatomicForceField(diatomicPotentials));
        anumm.mdRun(0.2, 200);
    }
}

