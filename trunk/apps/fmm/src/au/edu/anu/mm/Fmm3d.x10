package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class implements the Fast Multipole Method for electrostatic
 * calculations in three dimensions.
 * <p>
 * The 3D simulation space is divided in an octree of <code>numLevels</code> levels.
 * </p> 
 * @author milthorpe
 */
public class Fmm3d {
    /** The number of levels in the octree. */
    public val numLevels : Int;

    public var numBoxes : Int;
    
    /** The eight top-level boxes in the octree. */
    public val topLevelBoxes : ValRail[FmmParentBox];

    val atoms : Rail[Atom];

    public def this(numLevels : Int, atoms : Rail[Atom]) {
        this.numLevels = numLevels;
        this.topLevelBoxes = ValRail.make[FmmParentBox](8, (Int) => new FmmParentBox());
        this.atoms = atoms;
        setup();
    }
    
    public def calculateEnergy() : Double {
        var energy : Double = 0.0;

        /*
            TODO
            multipole expansion at the lowest level

            multipole expansion transfer from child to parent box.
            
            local expansion for each level.

            return far field potential energy & error check.
        */

        return energy;
    }

    def setup() {
        remainingLevels : Int = numLevels - 1;
        for ((i) in 0..7) {
            createChildBoxes(topLevelBoxes(i), remainingLevels);
        }
        Console.OUT.println("numLevels = " + numLevels + " numBoxes = " + numBoxes);

        iter : Iterator[Atom] = atoms.iterator();
        while (iter.hasNext()) {
            thisAtom : Atom = iter.next();
            // TODO assign atoms to boxes
        }
    }

    /** Recursively creates the octree of boxes. */
    def createChildBoxes(box : FmmParentBox, remainingLevels : Int) {
        if (remainingLevels == 0) {
            for ((i) in 0..7) {
                box.children.add(new FmmLeafBox(box));
                numBoxes++;
            }
        } else {
            for ((i) in 0..7) {
                lowerLevelBox : FmmParentBox = new FmmParentBox(box);
                box.children.add(lowerLevelBox);
                numBoxes++;
                createChildBoxes(lowerLevelBox, remainingLevels - 1);
            }   
        }     
    }
}

