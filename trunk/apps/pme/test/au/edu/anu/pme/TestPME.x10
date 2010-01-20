package au.edu.anu.pme;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * Tests the Distributed Particle Mesh Ewald implementation.
 * @author milthorpe
 */
public class TestPME {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);

    public static def main(args : Rail[String]!) {
        var numParticles : Int;
        var ewaldCoefficient : Double = 0.35;
        var gridSize : Int = 12;
        var splineOrder : Int = 4;
        if (args.length > 0) {
            numParticles = Int.parseInt(args(0));
            if (args.length > 1) {
                ewaldCoefficient = Double.parseDouble(args(1));
                if (args.length > 2) {
                    gridSize = Int.parseInt(args(2));
                    if (args.length > 3) {
                        splineOrder = Int.parseInt(args(3));
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPME numParticles [ewaldCoefficient] [gridSize] [splineOrder]");
            return;
        }

        if (splineOrder > gridSize) {
            Console.ERR.println("TestPME: splineOrder must not be greater than gridSize");
            return;
        }

        /* Assign particles to random locations within a -1..1 3D box, with unit charge (1/3 are negative). */
        val atoms : Rail[Atom] = ValRail.make[Atom](numParticles, (i : Int) => new Atom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val size = 2.0; // side length of cubic unit cell
        val edges = [new Vector3d(size, 0.0, 0.0), new Vector3d(0.0, size, 0.0), new Vector3d(0.0, 0.0, size)];
        val g = gridSize;
        val gridSizes = ValRail.make[Int](3, (Int) => g);
        val energy : Double = new PME(edges, gridSizes, atoms, splineOrder, ewaldCoefficient).calculateEnergy();
        Console.OUT.println("energy = " + energy);
    }

    static def randomUnit() : Double {
        return (R.nextDouble()) * 2.0;
    }
}

