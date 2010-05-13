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
        // distribute all atoms to all places
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        finish ateach (p1 in atoms) {
            val myAtoms = atoms(p1);
            for (p2 in atoms) {
                val otherAtoms = at(atoms.dist(p2)) {atoms(p2)};
                finish foreach ((j) in 0..otherAtoms.length-1) {
                    var myDirectEnergy : Double = 0.0;
                    val otherAtomCentre = at(atoms.dist(p2)) {otherAtoms(j).centre};
                    for ((i) in 0..myAtoms.length-1) {
                        val myAtom = myAtoms(i);
                        if (p1 != p2 || i != j) {
                            myDirectEnergy += myAtom.charge * otherAtoms(j).charge / otherAtomCentre.distance(myAtom.centre);
                        }
                    }
                    val myDirectEnergyFinal = myDirectEnergy;
                    // TODO this is slow because of lack of optimized atomic - XTENLANG-321
                    at(this) {atomic { directEnergy += myDirectEnergyFinal; }}
                }
            }
        }
       
        timer.stop(TIMER_INDEX_TOTAL);
        return directEnergy / 2.0;
    }
}
