package au.edu.anu.mm;

import x10.util.Random;
import x10x.vector.Point3d;

/**
 * Tests the FMM 3D implementation.
 * @author milthorpe
 */
public class TestFmm3d {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);

    public static def main(args: Rail[String]) {
        var numParticles : Int;
        var density : Double = 3.0;
        var numTerms : Int = 10;
        var wellSpaced : Int = 2;
        if (args.length > 0) {
            numParticles = Int.parseInt(args(0));
            if (args.length > 1) {
                density = Double.parseDouble(args(1));
                if (args.length > 2) {
                    numTerms = Int.parseInt(args(2));
                    if (args.length > 3) {
                        wellSpaced = Int.parseInt(args(3));
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestFmm3d numParticles [density] [numTerms] [wellSpaced]");
            return;
        }

        
        /* Assign particles to random locations within a -1..1 3D box, with unit charge (1/3 are negative). */
        val atoms : Rail[Atom] = ValRail.make[Atom](numParticles, (i : Int) => new Atom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val topLeftFront = new Point3d(-1.0, -1.0, -1.0);
        val energy : Double = new Fmm3d(density, numTerms, wellSpaced, topLeftFront, 2.0, atoms).calculateEnergy();
        Console.OUT.println("energy = " + energy);
    }

    static def randomUnit() : Double {
        return (R.nextDouble() - 0.5) * 2.0;
    }
}

