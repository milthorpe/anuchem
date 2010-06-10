package au.edu.anu.mm;

import x10.util.*;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * This subclass of Fmm3d extends the base FMM with periodic boundary
 * conditions.  In addition to the interactions within the unit cell, 
 * the  unit cell interacts with 3^k * 3^k * 3^k copies of itself in 
 * concentric shells of increasingly coarse-grained aggregate cells.
 * @see Lambert, Darden & Board (1996). "A Multipole-Based Algorithm 
 *  for Efficient Calculation of Forces and Potentials in Macroscopic 
 *  Periodic Assemblies of Particles". J Comp Phys 126 274-285
 * @author milthorpe
 */
public class PeriodicFmm3d extends Fmm3d {
    /** The number of concentric shells of copies of the unit cell. */
    public global val numShells : Int;

    public const TIMER_INDEX_MACROSCOPIC : Int = 7;
    /** 
     * A multi-timer for the several segments of a single getEnergy 
     * invocation, indexed by the constants above and in the superclass. 
     */
    public val timer = new Timer(8);

    /**
     * Initialises a periodic fast multipole method electrostatics 
     * calculation for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     * @param numShells the number of concentric shells of copies of the unit cell
     */
    public def this(density : Double, 
                    numTerms : Int,
                    topLeftFront : Point3d,
                    size : Double,  
                    numAtoms : Int,
                    atoms: DistArray[ValRail[MMAtom]](1),
                    numShells : Int) {
        // Periodic FMM always uses ws = 1
        // TODO is it possible to formulate for well-spaced > 1?
        super(density, numTerms, 1, topLeftFront, size, numAtoms, atoms);
        this.numShells = numShells;
    }

    public global def getTopLevel() : Int = 0;
    
    public def calculateEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);
        multipoleLowestLevel();
        // direct energy is independent of all subsequent steps of FMM
        val directEnergy = future {getDirectEnergy()};
        combineMultipoles();
        combineMacroscopicExpansions();
        transformToLocal();
        farFieldEnergy = getFarFieldEnergy();
        val totalEnergy = directEnergy() + farFieldEnergy;
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    /** 
     * Generate and includes macroscopic expansions, for concentric rings
     * of aggregates of copies of the unit cell.
     */
    def combineMacroscopicExpansions() {
        timer.start(TIMER_INDEX_MACROSCOPIC);
        // TODO distributed impl.
         val macroMultipoles = DistArray.make[MultipoleExpansion](Dist.makeBlock([0..Place.MAX_PLACES-1,0..numShells-1],0));
        finish ateach ((p) : Point in Dist.makeUnique(Place.places)) {
            //Console.OUT.println("shells at " + p + " top level = " + getTopLevel());
            val threeCube : Region(3) = [-1..1,-1..1,-1..1] as Region(3);
            val sevenCube : Region(3) = [-7..7,-7..7,-7..7] as Region(3);
            val topLevelBox = boxes(0)(0,0,0) as FmmBox!;
            // "shell 0" is the unit cell itself
            macroMultipoles(p, 0) = topLevelBox.multipoleExp;
            for (var shell: Int = 1; shell < numShells; shell++) {
                val sideLength = size * Math.pow(3.0, shell);
                //Console.OUT.println("shell = " + shell + ", sidelength = " + sideLength);
                macroMultipoles(p, shell) = new MultipoleExpansion(numTerms);
                for ((i,j,k) in threeCube) {
                    val translationVector = Vector3d(i * sideLength,
                                                     j * sideLength,
                                                     k * sideLength);
                    val translation = MultipoleExpansion.getOlm(translationVector, numTerms);
                    macroMultipoles(p, shell).translateAndAddMultipole(translation, macroMultipoles(p, shell-1));
                }
                //Console.OUT.println("final for " + shell + " = " + macroMultipoles(p, shell));
            }

            // now transform and add macroscopic multipoles to local expansion for top level box
            for (var shell: Int = 0; shell < numShells; shell++) {
                val sideLength = size * Math.pow(3.0, shell);
                for ((i,j,k) in sevenCube) {
                    if (Math.abs(i) > 1 || Math.abs(j) > 1 || Math.abs(k) > 1) {
                        // inner 27 boxes done at a lower level
                        val translationVector = Vector3d(i * sideLength,
                                                         j * sideLength,
                                                         k * sideLength);
                        val transform = LocalExpansion.getMlm(translationVector, numTerms) as LocalExpansion!;
                        topLevelBox.localExp.transformAndAddToLocal(transform, macroMultipoles(p, shell) as MultipoleExpansion!);
                    }
                }
            }
            Console.OUT.println("final for topLevel = " + topLevelBox.localExp);
        }
        timer.stop(TIMER_INDEX_MACROSCOPIC);
    }

}

