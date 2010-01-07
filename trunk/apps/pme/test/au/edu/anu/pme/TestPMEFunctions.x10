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
        val atoms : Rail[Atom] = ValRail.make[Atom](5, (i : Int) => new Atom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val topLeftFront = new Point3d(-1.0, -1.0, -1.0);
        val pme = new PME(topLeftFront, 2.0, atoms);
        val v = new Vector3d(0.6, -0.4, 0.33);
        val r = pme.getScaledFractionalCoordinates(v);
        Console.OUT.println(r);
    }

    static def randomUnit() : Double {
        return (R.nextDouble() - 0.5) * 2.0;
    }
}

