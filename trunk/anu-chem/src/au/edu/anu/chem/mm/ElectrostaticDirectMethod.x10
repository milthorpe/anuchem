package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.util.Timer;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class ElectrostaticDirectMethod {
    // TODO enum - XTENLANG-1118
    public const TIMER_INDEX_TOTAL : Int = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(6);

	private val atoms : ValRail[MMAtom!];

    // TODO should be shared local to getEnergy() - XTENLANG-404
    private var directEnergy : Double = 0.0;

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell
     */
    public def this(atoms : ValRail[MMAtom!]) {
        this.atoms = atoms;

        Console.OUT.println("Direct for " + atoms.length + " particles.");
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        // TODO multiplace
        finish foreach ((i) in 0..atoms.length-1) {
            var myDirectEnergy : Double = 0.0;
            for (var j : Int = 0; j < i; j++) {
                val distance = new Vector3d(atoms(j).centre.sub(atoms(i).centre as Tuple3d)).length();
                myDirectEnergy += atoms(i).charge * atoms(j).charge / distance;
            }
            // TODO this is slow because of lack of optimized atomic - XTENLANG-321
            atomic { directEnergy += myDirectEnergy; }
        }
       
        timer.stop(TIMER_INDEX_TOTAL);
        return directEnergy;
    }
}
