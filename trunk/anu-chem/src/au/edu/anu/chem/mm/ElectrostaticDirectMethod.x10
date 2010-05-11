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

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : DistArray[ValRail[MMAtom]](1);

    // TODO should be shared local to getEnergy() - XTENLANG-404
    private var directEnergy : Double = 0.0;

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell
     */
    public def this(atoms : DistArray[ValRail[MMAtom]](1)) {
        this.atoms = atoms;
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        finish ateach (p1 in atoms) {
            val localAtoms = atoms(p1);
            foreach ((i) in 0..localAtoms.length-1) {
                var myDirectEnergy : Double = 0.0;
                for (var j : Int = 0; j < i; j++) {
                    myDirectEnergy += localAtoms(i).charge * localAtoms(j).charge / localAtoms(j).centre.distance(localAtoms(i).centre);
                }
                val myDirectEnergyFinal = myDirectEnergy;
                // TODO this is slow because of lack of optimized atomic - XTENLANG-321
                at(this) {atomic { directEnergy += myDirectEnergyFinal; }}
            }
        }
       
        timer.stop(TIMER_INDEX_TOTAL);
        return directEnergy;
    }
}
