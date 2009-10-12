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

    /** All boxes in the octree division of space. */
    val boxes : Array[Box[FmmBox!]]{rank==2};

    val atoms : ValRail[Atom];

    /** A cache of transformations from multipole to local at the same level. */
    val multipoleTransforms : Array[LocalExpansion]{rank==4};

    /** A cache of multipole translations between parent box centres and child box centres. */
    val multipoleTranslations : Array[MultipoleExpansion]{rank==4};

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
                    atoms : ValRail[Atom]) {
        val numLevels = Math.max(2, (Math.log(atoms.length / density) / Math.log(8.0) + 1.0 as Int));
        this.numLevels = numLevels;
        var nBox : Int = 0;
        for ((i) in 2..numLevels) {
            nBox += Math.pow(8,i) as Int;
        }
        this.dimLowestLevelBoxes = Math.pow2(numLevels);
        Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + nBox);

        this.numTerms = numTerms;
        this.ws = ws;

        this.topLeftFront = topLeftFront;
        this.size = size;

        this.atoms = atoms;

        var boxRegion : Region{rank==2} = [0..63, 2..2];
        for ((i) in 3..numLevels) {
            rNextLevel : Region{rank==2} = [0..(Math.pow(8,i) as Int)-1, i..i];
            boxRegion = boxRegion || rNextLevel;
        }
        Console.OUT.println("boxes: " + boxRegion);

        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[Box[FmmBox!]](boxRegion);
    
        var wellSpacedLimit : Region(4) = [2..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransformRegion : Region(4) = wellSpacedLimit - ([2..numLevels,-ws..ws,-ws..ws,-ws..ws] as Region);
        this.multipoleTransforms = Array.make[LocalExpansion](multipoleTransformRegion);
        this.multipoleTranslations = Array.make[MultipoleExpansion]([2..numLevels, 0..1, 0..1, 0..1]);
    }
    
    public def calculateEnergy() : Double {

        // precompute multipole translations
        for (val(level,i,j,k) in multipoleTranslations.region) {
            dim : Int = Math.pow2(level);
            sideLength : Double = size / dim;
            translationVector : Vector3d = new Vector3d((i*2-1) * 0.5 * sideLength,
                                                             (j*2-1) * 0.5 * sideLength,
                                                             (k*2-1) * 0.5 * sideLength);
            multipoleTranslations(level, i, j, k) = MultipoleExpansion.getOlm(translationVector, numTerms);
        }

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
            val atom = atoms(i);
            val boxLocation = getLowestLevelBoxLocation(atom);
            val boxIndex = FmmBox.getBoxIndex(boxLocation, numLevels);
            val parentBox = getParentBox(boxIndex, numLevels);
            var box : Box[FmmBox!] = boxes(boxIndex, numLevels);
            if (box == null) {
                box = new Box[FmmBox!](new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox));
                boxes(boxIndex, numLevels) = box;
            }
            val leafBox = box.value as FmmLeafBox!;
            leafBox.atoms.add(atom);
            //Console.OUT.println("atoms(" + i + ") => box(" + boxIndex + ")");
            val boxCentre = leafBox.getCentre(size).sub(atom.centre);
            leafBox.multipoleExp.add(MultipoleExpansion.getOlm(atom.charge, boxCentre, numTerms));
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            for ((childIndex) in 0..(Math.pow(8,level) as Int)-1) {
                val bChild = boxes(childIndex,level);
                if (bChild != null) {
                    val child = bChild.value;
                    val bParent = child.parent;
                    val parentExp = at(bParent.value) {bParent.value.multipoleExp};
                    val shift = multipoleTranslations(level, (child.gridLoc(0)+1)%2, (child.gridLoc(1)+1)%2, (child.gridLoc(2)+1)%2);
                    parentExp.translateAndAddMultipole(shift, child.multipoleExp);
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
        // precompute child translations
        for (val(level,i,j,k) in multipoleTransforms.region) {
            dim : Int = Math.pow2(level);
            sideLength : Double = size / dim;
            translationVector : Vector3d = new Vector3d(i * sideLength,
                                                        j * sideLength,
                                                        k * sideLength);
            multipoleTransforms(level, i, j, k) = LocalExpansion.getMlm(translationVector, numTerms);
        } 

        var wellSep : Int = 0;
        var nearField : Int = 0;
        // top level (==2)
        for ((boxIndex1) in 0..63) {
            val bBox1 = boxes(boxIndex1, 2);
            if (bBox1 != null) {
                val box1 = bBox1.value;
                //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + 2 + ")");
                for ((boxIndex2) in 0..boxIndex1-1) { 
                    val bBox2 = boxes(boxIndex2, 2);
                    if (bBox2 != null) {
                        val box2 = bBox2.value;
                        //Console.OUT.println("... and box(" + boxIndex2 + "," + 2 + ")");
                        if (box2.wellSeparated(ws, box1)) {
                            val translation = getTranslationIndex(box1, box2);
                            val transform12 = multipoleTransforms(2, translation(0), translation(1), translation(2));
                            box2.localExp.transformAndAddToLocal(transform12, box1.multipoleExp);
                            val transform21 = multipoleTransforms(2, -translation(0), -translation(1), -translation(2));
                            box1.localExp.transformAndAddToLocal(transform21, box2.multipoleExp);
                            wellSep++;
                        } else if (numLevels==2) {
                            nearField++;
                        }
                    }
                }
            }
        }

        for ((level) in 3..numLevels) {
            for ((boxIndex1) in 0..(Math.pow(8,level) as Int)-1) {
                val bBox1 = boxes(boxIndex1, level);
                if (bBox1 != null) {
                    val box1 = bBox1.value;
                    val box1Parent = at(box1.parent.location) {box1.parent.value} as FmmBox!;
                    //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ")");
                    for ((boxIndex2) in 0..boxIndex1-1) { 
                        val bBox2 = boxes(boxIndex2, level);
                        if (bBox2 != null) {
                            val box2 = bBox2.value;
                            //Console.OUT.println("... and box(" + boxIndex2 + "," + level + ")");
                            if (!box2.parent.value.wellSeparated(ws, box1.parent.value)) {
                                if (box2.wellSeparated(ws, box1)) {
                                    val translation = getTranslationIndex(box1, box2);
                                    val transform12 = multipoleTransforms(level, translation(0), translation(1), translation(2));
                                    box2.localExp.transformAndAddToLocal(transform12, box1.multipoleExp);
                                    val transform21 = multipoleTransforms(level, -translation(0), -translation(1), -translation(2));
                                    box1.localExp.transformAndAddToLocal(transform21, box2.multipoleExp);
                                    wellSep++;
                                    
                                } else if (level==numLevels) {
                                    nearField++;
                                }
                            }
                        }
                    }
                    val shift = multipoleTranslations(level, box1.gridLoc(0)%2, box1.gridLoc(1)%2, box1.gridLoc(2)%2);
                    box1.localExp.translateAndAddLocal(shift, box1Parent.localExp);
                }
            }
        }
        Console.OUT.println("wellSep = " + wellSep + " nearField = " + nearField);
    }

    def getEnergy() : Double {
        // TODO n^2 calculation - to check - remove this
        /*
        var directEnergy : Double = 0.0;
        for ((i) in 0..(atoms.length - 1)) {
            for ((j) in 0..(i - 1)) {
                val pairEnergy : Double = atoms(j).pairEnergy(atoms(i));
                directEnergy += 2 * pairEnergy;
            }
        }
        Console.OUT.println("directEnergy = " + directEnergy);
        */

        var fmmEnergy : Double = 0.0;
        for ((boxIndex1) in 0..(Math.pow(8,numLevels) as Int)-1) { 
            val bBox1 = boxes(boxIndex1, numLevels);
            if (bBox1 != null) {
                val box1 = bBox1.value as FmmLeafBox!;
                for ((atomIndex1) in 0..box1.atoms.length()-1) {
                    val atom1 = box1.atoms(atomIndex1);
                    val box1Centre = atom1.centre.sub(box1.getCentre(size));
                    val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                    fmmEnergy += farFieldEnergy;

                    // direct calculation with all atoms in same box
                    for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                        val sameBoxAtom = box1.atoms(sameBoxAtomIndex) as Atom!;
                        val pairEnergy : Double = sameBoxAtom.pairEnergy(atom1);
                        fmmEnergy += 2 * pairEnergy;
                    }

                    // direct calculation with all atoms in non-well-separated boxes
                    for ((boxIndex2) in 0..boxIndex1-1) {
                        val bBox2 = boxes(boxIndex2, numLevels);
                        if (bBox2 != null) {
                            val box2 = bBox2.value as FmmLeafBox!;
                            if (!box1.wellSeparated(ws, box2)) {
                                for ((atomIndex2) in 0..box2.atoms.length()-1) {
                                    val atom2 = box2.atoms(atomIndex2) as Atom!;
                                    val pairEnergy : Double = atom2.pairEnergy(atom1);
                                    fmmEnergy += 2 * pairEnergy;
                                }
                            }
                        }
                    }
                }
            }
        }

        return fmmEnergy;
    }

    private def getLowestLevelBoxLocation(atom : Atom) : ValRail[Int]{length==3} {
        index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        return index;
    }

    private def getParentBox(childIndex : Int, childLevel : Int) : Box[FmmBox] {
        //Console.OUT.println("getParentBox(" + childIndex + ", " + childLevel + ")");
        if (childLevel == 2)
            return null;
        val parentIndex = getParentIndex(childIndex, childLevel);
        val parentLevel = childLevel - 1;
        var parent : Box[FmmBox] = boxes(parentIndex, parentLevel);
        if (parent == null) {
            parent = createNewParent(parentIndex, parentLevel);
        }
        return parent;
    }

    private def createNewParent(parentIndex : int, parentLevel : Int) : FmmBox {
        val grandparent = getParentBox(parentIndex, parentLevel);
        val parentLocation = getBoxLocation(parentIndex, parentLevel);
        val newParent = new FmmBox(parentLevel, parentLocation, numTerms, grandparent);
        boxes(parentIndex, parentLevel) = newParent;
        return newParent;
    }

    private def getBoxLocation(index : Int, level : Int) : ValRail[Int]{length==3} {
        dim : Int = Math.pow2(level) as Int;
        gridLoc : ValRail[Int]{length==3} = [ index / (dim * dim), (index / dim) % dim, index % dim ];
        return gridLoc;
    }

    private def getParentIndex(childIndex : Int, childLevel : Int) : Int {
        parentDim : Int = Math.pow2(childLevel-1) as Int;
        gridLoc : ValRail[Int] = getBoxLocation(childIndex, childLevel);
        return gridLoc(0) / 2 * parentDim * parentDim + gridLoc(1) / 2 * parentDim + gridLoc(2) / 2;
    }

    private def getTranslationIndex(box1 : FmmBox, box2 : FmmBox) : ValRail[Int]{length==3} {
        return [box1.gridLoc(0) - box2.gridLoc(0), box1.gridLoc(1) - box2.gridLoc(1), box1.gridLoc(2) - box2.gridLoc(2)];
    }
}

