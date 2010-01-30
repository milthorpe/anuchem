package au.edu.anu.pme;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

/**
 * Tests the Distributed Particle Mesh Ewald implementation - basic functions.
 * @author milthorpe
 */
public class TestPMEFunctions {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);

    public static def main(args : Rail[String]!) {
        val atoms = ValRail.make[MMAtom](3, (i : Int) => new MMAtom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val size = 2.0; // side length of cubic unit cell
        val edges = [new Vector3d(size, 0.0, 0.0), new Vector3d(0.0, size, 0.0), new Vector3d(0.0, 0.0, size)];
        val gridSize = ValRail.make[Int](3, (Int) => 12);
        val pme = new PME(edges, gridSize, atoms, 4, 0.35, 9.0);
        val v = new Vector3d(0.6, -0.4, 0.33);
        val r = pme.getScaledFractionalCoordinates(v);
        Console.OUT.println(r);

        val m404 = pme.bSpline(4, 0.4);
        Console.OUT.println("M_4(0.4) = " + m404);

        for (var a : int=0; a<atoms.length; a++) {
            Console.OUT.println(atoms(a));
        }

        val q = pme.getGriddedCharges();
        Console.OUT.println(q);
        var nonZero : Int = 0;
        for (p in q) {
            if (q(p) != Complex.ZERO) {
                nonZero++;
                //Console.OUT.println(p + " = " + q(p));
            }
        }
        Console.OUT.println("nonZero = " + nonZero);

        val B = pme.getBArray();
        var sum : Double = 0.0;
        for (p in B) {
            sum += B(p);
        }
        Console.OUT.println("sum B = " + sum);

        val C = pme.getCArray();
        for (p in C) {
            if (C(p) != 0.0) {
                Console.OUT.println(p + " = " + C(p));
            }
        }
    }

    static def randomUnit() : Double {
        return (R.nextDouble()) * 2.0;
    }
}

