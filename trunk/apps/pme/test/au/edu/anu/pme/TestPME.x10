package au.edu.anu.pme;

import x10.util.Random;
import x10.util.GrowableRail;
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

        val edges = [Vector3d(size, 0.0, 0.0), Vector3d(0.0, size, 0.0), Vector3d(0.0, 0.0, size)];
        val g = gridSize;
        val gridSizes = ValRail.make[Int](3, (Int) => g);

        Console.OUT.println("Testing PME for " + numParticles + " particles."
            + "\nBox edges: " + edges
            + "\nGrid size: " + gridSize
            + "\nspline order: " + splineOrder + " Beta: " + ewaldCoefficient + " Cutoff: " + cutoff);

        val atoms = generateAtoms(numParticles);
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

    /**
     * Generate an array of ValRails of MMAtoms, one ValRail for each
     * place.  PME assumes that the atoms have already been distributed. 
     */
    public static def generateAtoms(numAtoms : Int) : DistArray[ValRail[MMAtom]](1) {
        val tempAtoms = DistArray.make[GrowableRail[MMAtom]](Dist.makeUnique(Place.places), (Point) => new GrowableRail[MMAtom]());
        /* Assign particles to random locations within a small cubic area around the center of the simulation space, with unit charge (1/2 are negative). */
        finish for (var i : Int = 0; i < numAtoms; i++) {
            val x = randomUnit();
            val y = randomUnit();
            val z = randomUnit();
            val charge = i%2==0?1:-1;
            val p = getPlaceId(x, y, z);
            async (Place.places(p)) {
                val atom = new MMAtom(Point3d(x, y, z), 0.0, charge);
                atomic { (tempAtoms(p) as GrowableRail[MMAtom]!).add(atom); }
            }
        }
        val atoms = DistArray.make[ValRail[MMAtom]](Dist.makeUnique(Place.places), ((p) : Point) => (tempAtoms(p) as GrowableRail[MMAtom]!).toValRail());
        return atoms;
    }

    /** 
     * Gets the place ID to which to assign the given atom coordinates.
     * Currently just splits them up into slices by X coordinate.
     */
    public static safe def getPlaceId(x : Double, y : Double, z : Double) : Int {
        return ((x / size) * Place.MAX_PLACES) as Int;
    }

    static def randomUnit() : Double {
        val dub = at(R){R.nextDouble()};
        // uncomment to put a huge empty border around the particles
        //return dub * 2.0 - 1.0 + size / 2.0;
        return dub * size;
    }
}

