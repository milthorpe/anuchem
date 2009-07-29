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

    /** All boxes in the octree division of space. */
    val boxes : Array[FmmBox]{rank==2};

    val atoms : Rail[Atom];

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
                    atoms : Rail[Atom]) {
        this.numLevels = Math.max(2, (Math.log(atoms.length / density) / Math.log(8.0) + 1.0 as Int));
        var nBox : Int = 1;
        for ((i) in 2..numLevels) {
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

        var boxRegion : Region{rank==2} = [0..63, 2..2];
        for ((i) in 3..numLevels) {
            rNextLevel : Region{rank==2} = [0..(Math.pow(8,i) as Int)-1, i..i];
            boxRegion = boxRegion || rNextLevel;
        }
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[FmmBox](boxRegion);
        
        Console.OUT.println("boxes: " + boxes.region);
    
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
            atom : Atom = atoms(i);
            boxLocation : ValRail[Int]{length==3} = getLowestLevelBoxLocation(atom);
            boxIndex : Int = FmmBox.getBoxIndex(boxLocation, numLevels);
            parentBox : FmmParentBox = getParentBox(boxIndex, numLevels);
            var box : FmmLeafBox = boxes(boxIndex, numLevels) as FmmLeafBox;
            if (box == null) {
                box = new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox);
                if (parentBox != null) {
                    parentBox.children.add(box);
                }
                boxes(boxIndex, numLevels) = box;
            }
            box.atoms.add(atom);
            //Console.OUT.println("atoms(" + i + ") => box(" + boxIndex + ")");
            v : Tuple3d = box.getCentre(size).sub(atom.centre);
            olm : MultipoleExpansion = MultipoleExpansion.getOlm(atom.charge, v, numTerms);
            box.multipoleExp.add(olm);
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            for ((childIndex) in 0..(Math.pow(8,level) as Int)-1) {
                child : FmmBox = boxes(childIndex,level);
                if (child != null) {
                    val parent : FmmBox = child.parent;
                    val shift : MultipoleExpansion = multipoleTranslations(level, (child.location(0)+1)%2, (child.location(1)+1)%2, (child.location(2)+1)%2);
                    MultipoleExpansion.translateAndAddMultipole(shift, child.multipoleExp, parent.multipoleExp);
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
            box1 : FmmBox = boxes(boxIndex1, 2);
            if (box1 != null) {
                //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + 2 + ")");
                for ((boxIndex2) in 0..boxIndex1-1) { 
                    box2 : FmmBox = boxes(boxIndex2, 2);
                    if (box2 != null) {
                        //Console.OUT.println("... and box(" + 2 + "," + boxIndex2 + ")");
                        if (box1.wellSeparated(ws, box2)) {
                            translation : ValRail[Int]{length==3} = getTranslationIndex(box1, box2);
                                    transform12 : LocalExpansion = multipoleTransforms(2, translation(0), translation(1), translation(2));
                                    MultipoleExpansion.transformAndAddToLocal(transform12, box1.multipoleExp, box2.localExp);
                                    transform21 : LocalExpansion = multipoleTransforms(2, -translation(0), -translation(1), -translation(2));
                                    MultipoleExpansion.transformAndAddToLocal(transform21, box2.multipoleExp, box1.localExp);
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
                box1 : FmmBox = boxes(boxIndex1, level);
                if (box1 != null) {
                    //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ")");
                    for ((boxIndex2) in 0..boxIndex1-1) { 
                        box2 : FmmBox = boxes(boxIndex2, level);
                        if (box2 != null) {
                            //Console.OUT.println("... and box(" + level + "," + boxIndex2 + ")");
                            if (!box1.parent.wellSeparated(ws, box2.parent)) {
                                if (box1.wellSeparated(ws, box2)) {
                                    translation : ValRail[Int]{length==3} = getTranslationIndex(box1, box2);
                                    transform12 : LocalExpansion = multipoleTransforms(level, translation(0), translation(1), translation(2));
                                    MultipoleExpansion.transformAndAddToLocal(transform12, box1.multipoleExp, box2.localExp);
                                    transform21 : LocalExpansion = multipoleTransforms(level, -translation(0), -translation(1), -translation(2));
                                    MultipoleExpansion.transformAndAddToLocal(transform21, box2.multipoleExp, box1.localExp);
                                    wellSep++;
                                    
                                } else if (level==numLevels) {
                                    nearField++;
                                }
                            }
                        }
                    }
                    
                    shift : MultipoleExpansion = multipoleTranslations(level, box1.location(0)%2, box1.location(1)%2, box1.location(2)%2);
                    LocalExpansion.translateAndAddLocal(shift, box1.parent.localExp, box1.localExp);
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
                val pairEnergy : Double = pairEnergy(atoms(i), atoms(j));
                directEnergy += 2 * pairEnergy;
            }
        }
        Console.OUT.println("directEnergy = " + directEnergy);
        */

        var fmmEnergy : Double = 0.0;
        for ((boxIndex1) in 0..(Math.pow(8,numLevels) as Int)-1) { 
            box : FmmLeafBox = boxes(boxIndex1, numLevels) as FmmLeafBox;
            if (box != null) {
                for ((atomIndex1) in 0..box.atoms.length()-1) {
                    atom1 : Atom = box.atoms(atomIndex1);
                    v : Tuple3d = atom1.centre.sub(box.getCentre(size));
                    farFieldEnergy : Double = LocalExpansion.getPotential(atom1.charge, v, box.localExp);
                    fmmEnergy += farFieldEnergy;

                    // direct calculation with all atoms in same box
                    for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                        sameBoxAtom : Atom = box.atoms(sameBoxAtomIndex);
                        val pairEnergy : Double = pairEnergy(atom1, sameBoxAtom);
                        fmmEnergy += 2 * pairEnergy;
                    }

                    // direct calculation with all atoms in non-well-separated boxes
                    for ((boxIndex2) in 0..boxIndex1-1) {
                        box2 : FmmLeafBox = boxes(boxIndex2, numLevels) as FmmLeafBox;
                        if (box2 != null) {
                            if (!box.wellSeparated(ws, box2)) {
                                for ((atomIndex2) in 0..box2.atoms.length()-1) {
                                    atom2 : Atom = box2.atoms(atomIndex2);
                                    val pairEnergy : Double = pairEnergy(atom1, atom2);
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

    def pairEnergy(atom1 : Atom, atom2 : Atom) : Double {
        return atom1.charge * atom2.charge / (atom1.centre.distance(atom2.centre));
    }

    def getLowestLevelBoxLocation(atom : Atom) : ValRail[Int]{length==3} {
        index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        return index;
    }

    def getParentBox(childIndex : Int, childLevel : Int) : FmmParentBox {
        if (childLevel == 2)
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

    def getParentIndex(childIndex : Int, childLevel : Int) : Int {
        parentDim : Int = Math.pow2(childLevel-1) as Int;
        location : ValRail[Int] = getBoxLocation(childIndex, childLevel);
        return location(0) / 2 * parentDim * parentDim + location(1) / 2 * parentDim + location(2) / 2;
    }

    def getTranslationIndex(box1 : FmmBox, box2 : FmmBox) : ValRail[Int]{length==3} {
        return [box1.location(0) - box2.location(0), box1.location(1) - box2.location(1), box1.location(2) - box2.location(2)];
    }
}

