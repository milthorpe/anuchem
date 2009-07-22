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

        var boxRegion : Region{rank==2} = [0..8, 1..1];
        for ((i) in 2..numLevels) {
            rNextLevel : Region{rank==2} = [0..(Math.pow(8,i) as Int), i..i];
            boxRegion = boxRegion || rNextLevel;
        }
        
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[FmmBox](boxRegion);

        
        Console.OUT.println("boxes: " + boxes.region);
        
    }
    
    public def calculateEnergy() : Double {
        multipoleLowestLevel();

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
        combineMultipoles();
        /*
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                Console.OUT.println("boxes(" + i + "," + j + ") multipole = " + boxes(i,j).multipoleExp);
            }
        }
        */

        transformToLocal();
        /*
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                Console.OUT.println("boxes(" + i + "," + j + ") local = " + boxes(i,j).localExp);
            }
        }
        */

        return getEnergy();
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        for ((i) in 0..atoms.length-1) {
            atom : Atom = atoms(i);
            boxLocation : ValRail[Int]{length==3} = getLowestLevelBoxLocation(atom);
            boxIndex : Int = FmmBox.getBoxIndex(boxLocation, numLevels);
            parentBox : FmmParentBox = getParentBox(boxIndex, numLevels);
            var box : FmmLeafBox = boxes(boxIndex, numLevels) as FmmLeafBox;
            if (box == null) {
                box = new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox);
                parentBox.children.add(box);
                boxes(boxIndex, numLevels) = box;
            }
            v : Tuple3d = parentBox.getCentre(size).sub(atom.centre);
            olm : MultipoleExpansion = MultipoleExpansion.getOlm(atom.charge, v, numTerms);
            
            box.atoms.add(atom);
            //Console.OUT.println("atoms(" + i + ") => box(" + boxIndex + ")");
            box.multipoleExp.add(olm);
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 1; level--) {
            for ((childIndex) in 0..(Math.pow(8,level) as Int)) {
                child : FmmBox = boxes(childIndex,level);
                if (child != null) {
                    childCentre : Point3d = child.getCentre(size);
                    var parent : FmmBox = child.parent;
                    //Console.OUT.println("childIndex = " + childIndex + " parentIndex = " + parent.index());
                    v : Tuple3d = parent.getCentre(size).sub(childCentre);
                    MultipoleExpansion.translateAndAddMultipole(v, child.multipoleExp, parent.multipoleExp);
                    //Console.OUT.println(parent.multipoleExp);
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
                box1 : FmmBox = boxes(boxIndex1, level);
                if (box1 != null) {
                    //Console.OUT.println("transformToLocal: box(" + boxIndex1 + ")");
                    for ((boxIndex2) in 0..boxIndex1-1) { 
                        box2 : FmmBox = boxes(boxIndex2, level);
                        if (box2 != null) {
                            //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ") and box(" + level + "," + boxIndex2 + ")");
                            if (!box1.parent.wellSeparated(ws, box2.parent)) {
                                //Console.OUT.println("parents not well sep");
                                if (box1.wellSeparated(ws, box2)) {
                                    //Console.OUT.println("boxes well sep");
                                    wellSep++;
                                    boxCentre1 : Point3d = box1.getCentre(size);
                                    boxCentre2 : Point3d = box2.getCentre(size);
                                    v : Tuple3d = boxCentre2.sub(boxCentre1);
                                    MultipoleExpansion.transformAndAddToLocal(v, box1.multipoleExp, box2.localExp);
                                    MultipoleExpansion.transformAndAddToLocal(v.negate(), box2.multipoleExp, box1.localExp);
                                } else if (level==numLevels) {
                                    nearField++;
                                }
                            }
                        }
                    }
                    if (level > 2) {
                        v : Tuple3d = box1.getCentre(size).sub(box1.parent.getCentre(size));
                        LocalExpansion.translateAndAddLocal(v, box1.parent.localExp, box1.localExp);
                        //Console.OUT.println("after add parent: " + box1.localExp);
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
            var atomEnergy : Double = 0.0;
            for ((j) in 0..(atoms.length - 1)) {
                if (i != j) {
                    val pairEnergy : Double = pairEnergy(atoms(i), atoms(j));
                    atomEnergy += pairEnergy;
                    directEnergy += pairEnergy;
                }
            }
            Console.OUT.println("atomEnergy = " + atomEnergy);
        }

        Console.OUT.println("directEnergy = " + directEnergy);

        for ((boxIndex1) in 0..(Math.pow(8,numLevels) as Int)) { 
            box : FmmLeafBox = boxes(boxIndex1, numLevels) as FmmLeafBox;
            if (box != null) {
                for ((atomIndex1) in 0..box.atoms.length()-1) {
                    var atomEnergy : Double = 0.0;
                    atom1 : Atom = box.atoms(atomIndex1);
                    v : Tuple3d = atom1.centre.sub(box.getCentre(size));
                    //Console.OUT.println("atom(" + atomIndex1 + ") box(" + box.index() + ") v = " + v);
                    //Console.OUT.println("atom1 centre = " + atom1.centre + " box centre = " + box.getCentre(size));
                    //Console.OUT.println("localExp = " + box.localExp);
                    farFieldEnergy : Double = LocalExpansion.getPotential(atom1.charge, v, box.localExp);
                    //Console.OUT.println("farFieldEnergy = " + farFieldEnergy);
                    fmmEnergy += farFieldEnergy;
                    atomEnergy += farFieldEnergy;

                    for ((boxIndex2) in 0..boxIndex1-1) {
                        box2 : FmmLeafBox = boxes(boxIndex2, numLevels) as FmmLeafBox;
                        if (box2 != null) {
                            if (!box.wellSeparated(ws, box2)) {
                                //Console.OUT.println("box(" + boxIndex1 + ") and box(" + boxIndex2 + ") not well sep");
                                for ((atomIndex2) in 0..box2.atoms.length()-1) {
                                    atom2 : Atom = box2.atoms(atomIndex2);
                                    val pairEnergy : Double = pairEnergy(atom1, atom2);
                                    atomEnergy += pairEnergy;
                                    //Console.OUT.println("pairEnergy = " + pairEnergy);
                                    fmmEnergy += 2 * pairEnergy;
                                }
                            } else {
                                //Console.OUT.println("box(" + boxIndex1 + ") and box(" + boxIndex2 + ") well sep");
                            }
                        }
                    }
                    Console.OUT.println("atomEnergy = " + atomEnergy);
                }
            }
        }

        return fmmEnergy;
    }

    def pairEnergy(atom1 : Atom, atom2 : Atom) : Double {
        return atom1.charge * atom2.charge / (atom1.centre.distance(atom2.centre));
    }

    def getLowestLevelBoxLocation(atom : Atom) : ValRail[Int]{length==3} {
        //Console.OUT.println(atom.centre);
        index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        //Console.OUT.println(index(0) + " " + index(1) + " " + index(2));
        return index;
    }

    def getParentBox(childIndex : Int, childLevel : Int) : FmmParentBox {
        //Console.OUT.println("getParentBox(" + childIndex + ", " + childLevel);
        if (childLevel == 1)
            return null;
        parentIndex : Int = getParentIndex(childIndex, childLevel);
        parentLevel : Int = childLevel - 1;
        parentLocation : ValRail[Int]{length==3} = getBoxLocation(parentIndex, parentLevel);
        var parent : FmmParentBox = boxes(parentIndex, parentLevel) as FmmParentBox;
        if (parent == null) {
            grandparent : FmmParentBox = getParentBox(parentIndex, parentLevel);
            parent = new FmmParentBox(parentLevel, parentLocation, numTerms, grandparent);
            if (grandparent != null)
                grandparent.children.add(parent);
            boxes(parentIndex, parentLevel) = parent;
        }
        return parent;
    }

    def getBoxLocation(index : Int, level : Int) : ValRail[Int]{length==3} {
        dim : Int = Math.pow2(level) as Int;
        location : ValRail[Int]{length==3} = [ index / (dim * dim), (index / dim) % dim, index % dim ];
        return location;
    }

    /**
     * Returns true if <code>boxes(boxIndex1)</code> and 
     * <code>boxes(boxIndex2)</code> are well separated i.e. whether
     * there are at least <code>ws</code> boxes separating them.
     */
    def wellSeparated(ws : Int, boxIndex1 : Int, boxIndex2 : Int, level : Int) : Boolean {
        if (level < 2)
            return false;
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
        //Console.OUT.println("childLoc = " + location(0) + " " + location(1) + " " + location(2));
        return location(0) / 2 * parentDim * parentDim + location(1) / 2 * parentDim + location(2) / 2;
    }
}

