package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * Tests ANU-Mechanics with a simple Morse oscillator: HF.
 * A test problem borrowed from Chapter 1 of
 * Herman J.C. Berendsen, "Simulating the Physical World"
 * Cambridge University Press, 2007
 * @author milthorpe
 */
public class TestMorseOscillator {
    public static def main(args : Rail[String]!) {

        Console.OUT.println("Testing simple HF morse oscillator.");

        // equilibrium bond length is 0.09169nm; start with displacement of 0.3nm
        val hydrogen = new MMAtom("H", new Point3d(0.0, 0.0, 0.0), 1.0079, 0.0);
        val fluorine = new MMAtom("F", new Point3d(0.03, 0.0, 0.0), 18.9984, 0.0);
        val atoms = Rail.make[MMAtom](2);
        atoms(0) = hydrogen;
        atoms(1) = fluorine;
        val distAtoms = DistArray.make[ValRail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms as ValRail[MMAtom];

        val diatomicMolecules : ValRail[Pair[MMAtom,MMAtom]] = null;
        val anum = new Anum(distAtoms, new MorseForceField(diatomicMolecules));
        anum.mdRun(1.0, 200);
    }
}

