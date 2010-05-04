package au.edu.anu.pme;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.util.Timer;

/**
 * Tests the Distributed Particle Mesh Ewald implementation.
 * @author milthorpe
 */
public class TestPME {
    private const RANDOM_SEED = 10101110L;
    private static val R = new Random(RANDOM_SEED);
    /* side length of cubic unit cell in Angstroms */
    private const size = 80.0;

    public static def main(args : Rail[String]!) {
        var numParticles : Int;
        var ewaldCoefficient : Double = 0.35;
        var cutoff : Double = 10.0;
        var gridSize : Int = 64;
        var splineOrder : Int = 4;
        if (args.length > 0) {
            numParticles = Int.parseInt(args(0));
            if (args.length > 1) {
                ewaldCoefficient = Double.parseDouble(args(1));
                if (args.length > 2) {
                    cutoff = Double.parseDouble(args(2));
                    if (args.length > 3) {
                        gridSize = Int.parseInt(args(3));
                        if (args.length > 4) {
                            splineOrder = Int.parseInt(args(4));
                        }
                    }
                }
            }
        } else {
            Console.ERR.println("usage: TestPME numParticles [ewaldCoefficient] [cutoff] [gridSize] [splineOrder]");
            return;
        }

        if (splineOrder > gridSize) {
            Console.ERR.println("TestPME: splineOrder must not be greater than gridSize");
            return;
        }

        /* Assign particles to random locations within a small cubic area around the center of the simulation space, with unit charge (1/2 are negative). */
        val atoms = ValRail.make[MMAtom!](numParticles, (i : Int) => new MMAtom(Point3d(randomUnit(R) + size/2.0, randomUnit(R) + size/2.0, randomUnit(R) + size/2.0), 0.0, i%2==0?1:-1));
        val edges = [Vector3d(size, 0.0, 0.0), Vector3d(0.0, size, 0.0), Vector3d(0.0, 0.0, size)];
        val g = gridSize;
        val gridSizes = ValRail.make[Int](3, (Int) => g);
        val pme = new PME(edges, gridSizes, atoms, splineOrder, ewaldCoefficient, cutoff);
        val energy = pme.getEnergy();
        Console.OUT.println("energy = " + energy);

        logTime("Divide",            PME.TIMER_INDEX_DIVIDE,        pme.timer);
        logTime("Direct",            PME.TIMER_INDEX_DIRECT,        pme.timer);
        logTime("Self energy",       PME.TIMER_INDEX_SELF,          pme.timer);
        logTime("Grid charges",      PME.TIMER_INDEX_GRIDCHARGES,   pme.timer);
        logTime("Inverse FFT",       PME.TIMER_INDEX_INVFFT,        pme.timer);
        logTime("ThetaRecConvQ",     PME.TIMER_INDEX_THETARECCONVQ, pme.timer);
        logTime("Reciprocal energy", PME.TIMER_INDEX_RECIPROCAL,    pme.timer);
        logTime("Total",             PME.TIMER_INDEX_TOTAL,         pme.timer);

        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("cf. Direct calculation", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
        // direct error comparison is only useful if there is a huge empty border around the particles
        val error = directEnergy - energy;
        Console.OUT.println("direct = " + directEnergy + " error = " + error + " relative error = " + Math.abs(error) / Math.abs(energy));

    }

    private static def logTime(desc : String, timerIndex : Int, timer : Timer!) {
        Console.OUT.printf(desc + " (one cycle): %g seconds\n", (timer.total(timerIndex) as Double) / 1e9);
    }

    static def randomUnit(R : Random) : Double {
        val dub = at(R){R.nextDouble()};
        // uncomment to put a huge empty border around the particles
        //return (dub) * 2.0 - 1.0;
        return ((dub)-0.5) * size;
    }
}

