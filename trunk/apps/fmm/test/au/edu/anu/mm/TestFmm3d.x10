package au.edu.anu.mm;

import x10.util.Random;
import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed FMM 3D implementation.
 * @author milthorpe
 */
public class TestFmm3d {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);

    public static def main(args : Rail[String]!) {
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
        val atoms = ValRail.make[MMAtom!](numParticles, (i : Int) => new MMAtom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1));
        val fmm3d = new Fmm3d(density, numTerms, wellSpaced, 2.0, atoms);
        val timer = new Timer(1);
        timer.start(0);
        val energy = fmm3d.calculateEnergy();
        timer.stop(0);
        Console.OUT.println("energy = " + energy);
        Console.OUT.printf("Time (one cycle): %g seconds", (timer.total(0) as Double) / 1e9);
    }

    static def randomUnit() : Double {
        val dub = at(R){R.nextDouble()};
        return (dub - 0.5) * 2.0;
    }
}

