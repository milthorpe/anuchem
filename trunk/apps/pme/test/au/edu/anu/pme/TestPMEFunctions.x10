package au.edu.anu.pme;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * Tests the Distributed Particle Mesh Ewald implementation - basic functions.
 * @author milthorpe
 */
public class TestPMEFunctions {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);

    public static def main(args : Rail[String]!) {
        val atoms = ValRail.make[Atom](2, (i : Int) => new Atom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val gridSize = ValRail.make[Int](3, (Int) => 12);
        val pme = new PME(2.0, gridSize, atoms);
        val v = new Vector3d(0.6, -0.4, 0.33);
        val r = pme.getScaledFractionalCoordinates(v);
        Console.OUT.println(r);

        val m404 = pme.bSpline(4, 0.4);
        Console.OUT.println("M_4(0.4) = " + m404);

        for (var a : int=0; a<atoms.length; a++) {
            Console.OUT.println(atoms(a));
        }

        val q = pme.getGriddedCharges(2);
        Console.OUT.println(q);
        for (p in q) {
            if (q(p) != 0.0) {
                Console.OUT.println(p + " = " + q(p));
            }
        }

        val B = pme.getBArray(4);
        for (p in B) {
            if (B(p) != 0.0) {
                Console.OUT.println(p + " = " + B(p));
            }
        }
    }

    static def randomUnit() : Double {
        return (R.nextDouble()) * 2.0;
    }
}

