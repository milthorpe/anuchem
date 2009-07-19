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
        var numLevels : Int;
        var numParticles : Int;
        var numTerms : Int;
        if (args.length >= 1) {
            numParticles = Int.parseInt(args(0));
            numTerms = Int.parseInt(args(1));
            numLevels = Int.parseInt(args(2));
        } else {
            numParticles = 1000;
            numTerms = 5;
            numLevels = 2;
            //Console.ERR.println("usage: TestFmm3d <numParticles> <numTerms> <numLevels>");
        }

        
        /* Assign particles to random locations within a -1..1 3D box, with unit charge (1/3 are negative). */
        val atoms : Rail[Atom] = ValRail.make[Atom](numParticles, (i : Int) => new Atom("H", new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==0?1:-1));
        val topLeftFront = new Point3d(-1.0, -1.0, -1.0);
        val energy : Double = new Fmm3d(numLevels, numTerms, topLeftFront, 2.0, atoms).calculateEnergy();
    }

    static def randomUnit() : Double {
        return (R.nextDouble() - 0.5) * 2.0;
    }
}

