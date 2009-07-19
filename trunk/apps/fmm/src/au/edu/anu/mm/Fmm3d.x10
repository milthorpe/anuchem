package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;

/**
 * This class implements the Fast Multipole Method for electrostatic
 * calculations in a cubic simulation space.
 * <p>
 * The 3D simulation space is divided in an octree of <code>numLevels</code> levels.
 * </p> 
 * @author milthorpe
 */
public class Fmm3d {
    /** The number of levels in the octree. */
    public val numLevels : Int;

    /** The number of lowest level boxes along one side of the (cubic) 3D space. */
    public val dimLowestLevelBoxes : Int;

    /** The cartesian location of the top-left-front corner of the simulation cube. */
    public val topLeftFront : Point3d;

    /** The length of a side of the simulation cube. */
    public val size : Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public val numTerms : Int;

    /** The maximum number of boxes in the tree, if there are no empty boxes at the lowest level. */
    val maxBoxes : Int;

    /** The total number of non-empty boxes at the lowest level. */
    var numLeafBoxes : Int;

    /** All boxes in the octree division of space. */
    val boxes : Array[FmmBox]{rank==2};

    val atoms : Rail[Atom];
    val atomBoxIndices : Rail[Int];

    public def this(numLevels : Int, numTerms : Int, topLeftFront : Point3d, size : Double, atoms : Rail[Atom]) {
        this.numLevels = numLevels;
        var nBox : Int = 1;
        for ((i) in 1..numLevels) {
            nBox += Math.pow(8,i) as Int;
        }
        this.numTerms = numTerms;
        this.maxBoxes = nBox;   
        this.dimLowestLevelBoxes = Math.pow2(numLevels);
        this.topLeftFront = topLeftFront;
        this.size = size;
        this.atoms = atoms;
        this.atomBoxIndices = ValRail.make[Int](atoms.length, (i : Int) => getLowestLevelBoxIndex(atoms(i)));
        var boxRegion : Region{rank==2} = [0..8, 1..1];
        for ((i) in 2..numLevels) {
            rNextLevel : Region{rank==2} = [0..(Math.pow(8,i) as Int), i..i];
            boxRegion = boxRegion || rNextLevel;
        }
        
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[FmmBox](boxRegion);
    }
    
    public def calculateEnergy() : Double {
        var energy : Double = 0.0;

        Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + maxBoxes);
        Console.OUT.println("boxes: " + boxes.region);

        /*
        for ((i) in 0..atoms.length()-1) {
            Console.OUT.println("atom(" + i + ") index = " + atomBoxIndices(i));
        }
        */

        multipoleLowestLevel();

        combineMultipoles();

        var nonEmpty : Int = 0;
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                nonEmpty++;
                Console.OUT.println("boxes(" + i + "," + j + ") multipole = " + boxes(i,j).multipoleExp);
            }
        }
        Console.OUT.println("nonEmpty = " + nonEmpty);
        
        /*
            TODO
            
            local expansion for each level.

            return far field potential energy & error check.
        */

        return energy;
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        for ((i) in 0..atoms.length-1) {
            atom : Atom = atoms(i);
            boxIndex : Int = atomBoxIndices(i);
            v : Tuple3d = getBoxCentre(boxIndex, numLevels).sub(atom.centre);
            olm : MultipoleExpansion = MultipoleExpansion.getOlm(atom.charge, v, numTerms);
            var box : FmmBox = boxes(boxIndex, numLevels);
            if (box == null) {
                box = new FmmBox(numTerms);
                boxes(boxIndex, numLevels) = box;
            }
            box.multipoleExp.add(olm);
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level >= 2; level--) {
            for ((childIndex) in 0..(Math.pow(8,level) as Int)) {
                child : FmmBox = boxes(childIndex,level);
                if (child != null) {
                    childCentre : Point3d = getBoxCentre(childIndex, level);
                    parentIndex : Int = getParentIndex(childIndex, level);
                    //Console.OUT.println("childIndex = " + childIndex + " parentIndex = " + parentIndex);
                    var parent : FmmBox = boxes(parentIndex, level-1);
                    if (parent == null) {
                        //Console.OUT.println("assigning (" + parentIndex + "," + (level-1) + ")");
                        parent = new FmmBox(numTerms);
                        boxes(parentIndex, level-1) = parent;
                    }
                    parentCentre : Point3d = getBoxCentre(parentIndex, level-1);
                    //Console.OUT.println(parentCentre);
                    v : Tuple3d = parentCentre.sub(childCentre);
                    MultipoleExpansion.translateAndAddMultipole(v, child.multipoleExp, parent.multipoleExp);
                }
            }
        }
    }

    def getLowestLevelBoxLocation(atom : Atom) : Rail[Int] {
        //Console.OUT.println(atom.centre);
        index : ValRail[Int] = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        //Console.OUT.println(index(0) + " " + index(1) + " " + index(2));
        return index;
    }

    def getLowestLevelBoxIndex(atom : Atom) {
        index : ValRail[Int] = getLowestLevelBoxLocation(atom);
        return index(0) * dimLowestLevelBoxes * dimLowestLevelBoxes + index(1) * dimLowestLevelBoxes + index(2);
    }

    def getBoxLocation(index : Int, level : Int) : ValRail[int] {
        dim : Int = Math.pow2(level) as Int;
        location : ValRail[Int] = [ index / (dim * dim), (index / dim) % dim, index % dim ];
        return location;
    }

    def getParentIndex(childIndex : Int, childLevel : Int) : Int {
        parentDim : Int = Math.pow2(childLevel-1) as Int;
        location : ValRail[Int] = getBoxLocation(childIndex, childLevel);
        return location(0) / 2 * parentDim * parentDim + location(1) / 2 * parentDim + location(2) / 2;
    }
    
    def getBoxCentre(boxIndex : Int, level : Int) {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        return new Point3d( (boxIndex / (dim * dim) + 0.5) * sideLength - (0.5 * size),
                                (boxIndex % (dim * dim ) / (dim) + 0.5) * sideLength - (0.5 * size),
                                (boxIndex % dim + 0.5) * sideLength - (0.5 * size));
    }
}

