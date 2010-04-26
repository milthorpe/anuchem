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

        // equilibrium bond length is 0.09169nm; start with displacement of 0.1nm
        val hydrogen = new MMAtom("H", Point3d(0.1, 0.0, 0.0), 1.0079, 0.0);
        val fluorine = new MMAtom("F", Point3d(0.0, 0.0, 0.0), 18.9984, 0.0);
        val atoms = Rail.make[MMAtom](2);
        atoms(0) = hydrogen;
        atoms(1) = fluorine;
        val distAtoms = DistArray.make[ValRail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms as ValRail[MMAtom];

        val diatomicPotentials : ValRail[DiatomicPotential] = ValRail.make[DiatomicPotential](1, 
           (Int) => new DiatomicMorsePotential(hydrogen, fluorine, 0.09169, 569.87));
        val anum = new Anum(distAtoms, new DiatomicForceField(diatomicPotentials));
        anum.mdRun(0.2, 200);
    }
}

