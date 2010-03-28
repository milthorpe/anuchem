package au.edu.anu.mm;

import x10.util.Random;
import x10.util.GrowableRail;
import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * Tests the distributed FMM 3D implementation.
 * @author milthorpe
 */
public class TestFmm3d {
    private const RANDOM_SEED = 10101110L;
    private const size = 2.0;
    private static val R = new Random(RANDOM_SEED);
    private const X_SLICE = size / Place.MAX_PLACES;

    public static def main(args : Rail[String]!) {
        var numAtoms : Int;
        var density : Double = 60.0;
        var numTerms : Int = 10;
        var wellSpaced : Int = 2;
        if (args.length > 0) {
            numAtoms = Int.parseInt(args(0));
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
            Console.ERR.println("usage: TestFmm3d numAtoms [density] [numTerms] [wellSpaced]");
            return;
        }

        Console.OUT.println("Testing FMM for " + numAtoms 
                          + " atoms, target density = " + density
                          + " numTerms = " + numTerms
                          + " wellSpaced param = " + wellSpaced);
        


        val fmm3d = new Fmm3d(density, numTerms, wellSpaced, 2.0, numAtoms, generateAtoms(numAtoms));
        val energy = fmm3d.calculateEnergy();
        
        Console.OUT.println("energy = " + energy);

        val timer = fmm3d.timer;
        logTime("Multipole", Fmm3d.TIMER_INDEX_MULTIPOLE, fmm3d.timer);
        logTime("Direct",    Fmm3d.TIMER_INDEX_DIRECT,    fmm3d.timer);
        logTime("Combine",   Fmm3d.TIMER_INDEX_COMBINE,   fmm3d.timer);
        logTime("Transform", Fmm3d.TIMER_INDEX_TRANSFORM, fmm3d.timer);
        logTime("Far field", Fmm3d.TIMER_INDEX_FARFIELD,  fmm3d.timer);
        logTime("Total",     Fmm3d.TIMER_INDEX_TOTAL,     fmm3d.timer);
    }

    private static def logTime(desc : String, timerIndex : Int, timer : Timer!) {
        Console.OUT.printf(desc + " (one cycle): %g seconds\n", (timer.total(timerIndex) as Double) / 1e9);
    }

    /**
     * Generate an array of ValRails of MMAtoms, one ValRail for each
     * place.  FMM assumes that the atoms have already been distributed. 
     */
    public static def generateAtoms(numAtoms : Int) : Array[ValRail[MMAtom]](1) {
        val tempAtoms = Array.make[GrowableRail[MMAtom]](Dist.makeUnique(Place.places), (Point) => new GrowableRail[MMAtom]());
        /* Assign Atoms to random locations within a -1..1 3D box, with unit charge (1/3 are negative). */
        for (var i : Int = 0; i < numAtoms; i++) {
            val atom = new MMAtom(new Point3d(randomUnit(), randomUnit(), randomUnit()), i%3==4?1:-1);
            val p = getPlaceId(atom);
            at (Place.places(p)) {
                val remoteAtom = new MMAtom(atom);
                (tempAtoms(p) as GrowableRail[MMAtom]!).add(remoteAtom);
            }
        }
        val atoms = Array.make[ValRail[MMAtom]](Dist.makeUnique(Place.places), ((p) : Point) => (tempAtoms(p) as GrowableRail[MMAtom]!).toValRail());
        return atoms;
    }

    public static safe def getPlaceId(atom : MMAtom!) : Int {
        return ((atom.centre.i + 0.5 * size) / size * Place.MAX_PLACES) as Int;
    }

    static def randomUnit() : Double {
        val dub = at(R){R.nextDouble()};
        return dub * (size) - 0.5 * size;
    }
}

