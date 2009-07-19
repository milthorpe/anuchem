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

    /** The well-separatedness parameter ws. */
    public val ws : Int;

    /** The maximum number of boxes in the tree, if there are no empty boxes at the lowest level. */
    val maxBoxes : Int;

    /** The total number of non-empty boxes at the lowest level. */
    var numLeafBoxes : Int;

    /** All boxes in the octree division of space. */
    val boxes : Array[FmmBox]{rank==2};

    val atoms : Rail[Atom];
    val atomBoxIndices : Rail[Int];

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param topLeftFront cartesian location of the top-left-front corner of the simulation cube
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     */
    public def this(density : Double, 
                    numTerms : Int,
                    ws : Int,
                    topLeftFront : Point3d,
                    size : Double,
                    atoms : Rail[Atom]) {
        this.numLevels = Math.max(2, (Math.log(atoms.length / density) / Math.log(8.0) + 1.0 as Int));
        var nBox : Int = 1;
        for ((i) in 1..numLevels) {
            nBox += Math.pow(8,i) as Int;
        }
        this.maxBoxes = nBox;
        this.dimLowestLevelBoxes = Math.pow2(numLevels);
        Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + maxBoxes);

        this.numTerms = numTerms;
        this.ws = ws;

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

        /*
        Console.OUT.println("boxes: " + boxes.region);
        for ((i) in 0..atoms.length()-1) {
            Console.OUT.println("atom(" + i + ") index = " + atomBoxIndices(i));
        }
        */
    }
    
    public def calculateEnergy() : Double {
        multipoleLowestLevel();

        combineMultipoles();

        /*
        var nonEmpty : Int = 0;
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                nonEmpty++;
                Console.OUT.println("boxes(" + i + "," + j + ") multipole = " + boxes(i,j).multipoleExp);
            }
        }
        Console.OUT.println("nonEmpty = " + nonEmpty);
        */

        transformToLocal();
        
        return getEnergy();
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

    /**
     * Starting at the top level, for each box, transform multipole
     * expansions for all well-separated boxes (for which the parent
     * box is not also well-separated) to local expansions for this 
     * box.  Translate the local expansion for each parent down to all
     * non-empty child boxes.
     */
    def transformToLocal() {
        var wellSep : Int = 0;
        var nearField : Int = 0;
        for ((level) in 2..numLevels) {
            for ((boxIndex1) in 1..(Math.pow(8,level) as Int)) {
                for ((boxIndex2) in 0..boxIndex1-1) {
                    box1 : FmmBox = boxes(boxIndex1, level);
                    box2 : FmmBox = boxes(boxIndex2, level);
                    if (box1 != null && box2 != null) {
                        parentIndex1 : Int = getParentIndex(boxIndex1, level);
                        parentIndex2 : Int = getParentIndex(boxIndex2, level);
                        if (!wellSeparated(ws, parentIndex1, parentIndex2, level-1)) {
                            if (wellSeparated(ws, boxIndex1, boxIndex2, level)) {
                                //Console.OUT.println("boxes " + boxIndex1 + " and " + boxIndex2);
                                wellSep++;
                                boxCentre1 : Point3d = getBoxCentre(boxIndex1, level);
                                boxCentre2 : Point3d = getBoxCentre(boxIndex2, level);
                                v : Tuple3d = boxCentre1.sub(boxCentre2);
                                MultipoleExpansion.transformAndAddToLocal(v, box1.multipoleExp, box2.localExp);
                                MultipoleExpansion.transformAndAddToLocal(v.negate(), box2.multipoleExp, box1.localExp);
                            } else if (level==numLevels) {
                                nearField++;
                            }
                        }
                    }
                }
            }
        }
        Console.OUT.println("wellSep = " + wellSep + " nearField = " + nearField);
    }

    def getEnergy() : Double {
        var fmmEnergy : Double = 0.0;
        var directEnergy : Double = 0.0;

        // TODO n^2 calculation - to check - remove this
        for ((i) in 0..(atoms.length - 1)) {
            for ((j) in 0..(atoms.length - 1)) {
                if (i != j) {
                    val pairEnergy : Double = pairEnergy(atoms(i), atoms(j));
                    directEnergy += pairEnergy;
                    if (!wellSeparated(ws, atomBoxIndices(i), atomBoxIndices(j), numLevels)) {
                        // only add direct pair energy for particles in non-well-separated boxes
                        fmmEnergy += pairEnergy;
                    }
                }
            }
        }

        Console.OUT.println(directEnergy);

        for ((i) in 0..(atoms.length - 1)) {
            atom : Atom = atoms(i);
            box : FmmBox = boxes(atomBoxIndices(i), numLevels);
            v : Tuple3d = getBoxCentre(atomBoxIndices(i), numLevels).sub(atom.centre);
            farFieldEnergy : Double = LocalExpansion.getPotential(atom.charge, v, box.localExp);
            fmmEnergy += farFieldEnergy;
        }

        return fmmEnergy;
    }

    def pairEnergy(atom1 : Atom, atom2 : Atom) : Double {
        return atom1.charge * atom2.charge / (atom1.centre.distance(atom2.centre));
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

    /**
     * Returns true if <code>boxes(boxIndex1)</code> and 
     * <code>boxes(boxIndex2)</code> are well separated i.e. whether
     * there are at least <code>ws</code> boxes separating them.
     */
    def wellSeparated(ws : Int, boxIndex1 : Int, boxIndex2 : Int, level : Int) : Boolean {
        if (boxIndex1 == boxIndex2)
            return false;
        loc1 : ValRail[Int] = getBoxLocation(boxIndex1, level);
        loc2 : ValRail[Int] = getBoxLocation(boxIndex2, level);
        // TODO can do reduction on a Rail?
        return Math.abs(loc1(0) - loc2(0)) > ws 
            || Math.abs(loc1(1) - loc2(1)) > ws 
            || Math.abs(loc1(2) - loc2(2)) > ws;
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

