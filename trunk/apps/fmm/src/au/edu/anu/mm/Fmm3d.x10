package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;

/**
 * This class implements the Fast Multipole Method for electrostatic
 * calculations in a cubic simulation space centred at the origin.
 * <p>
 * The 3D simulation space is divided in an octree of <code>numLevels</code> levels.
 * </p> 
 * @author milthorpe
 */
public class Fmm3d {
    /** The number of levels in the octree. */
    public global val numLevels : Int;

    /** The number of lowest level boxes along one side of the cube. */
    public global val dimLowestLevelBoxes : Int;

    /** The side length of the cube. */
    public global val size : Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public global val numTerms : Int;

    /** The well-separatedness parameter ws. */
    public global val ws : Int;

    /** All boxes in the octree division of space. */
    global val boxes : Array[FmmBox](2);

    val atoms : ValRail[MMAtom!];

    /** A cache of transformations from multipole to local at the same level. */
    global val multipoleTransforms : Array[LocalExpansion](5);

    /** A cache of multipole translations between parent box centres and child box centres. */
    global val multipoleTranslations : Array[MultipoleExpansion](5);

    // TODO use shared local variable within getEnergy() - XTENLANG-404
    var fmmEnergy : Double = 0.0;

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     */
    public def this(density : Double, 
                    numTerms : Int,
                    ws : Int,
                    size : Double,
                    atoms : ValRail[MMAtom!]) {
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

        this.size = size;

        this.atoms = atoms;

        var boxRegion : Region(2) = [0..63, 2..2];
        var boxDistribution : Dist(2) = Dist.makeBlock(boxRegion, 0);
        for ((i) in 3..numLevels) {
            val nextLevelRegion : Region(2) = [0..(Math.pow(8,i) as Int)-1, i..i];
            val nextLevelDist = Dist.makeBlock(nextLevelRegion, 0);
            boxDistribution = boxDistribution || nextLevelDist;
        }
        Console.OUT.println("boxDistribution: " + boxDistribution);
        
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[FmmBox](boxDistribution);

        // two special arrays distributed to all places (this is done by replicating the first index in a cyclic dist across Place.PLACES)
        var wellSpacedLimit : Region(5) = [0..Place.MAX_PLACES-1,2..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransformRegion : Region(5) = wellSpacedLimit - ([0..Place.MAX_PLACES-1,2..numLevels,-ws..ws,-ws..ws,-ws..ws] as Region);
        //Console.OUT.println("multipoleTransformRegion = " + multipoleTransformRegion);
        this.multipoleTransforms = Array.make[LocalExpansion](Dist.makeCyclic(multipoleTransformRegion,0));
        if (numLevels >= 3) {
            this.multipoleTranslations = Array.make[MultipoleExpansion](Dist.makeCyclic([0..Place.MAX_PLACES-1,3..numLevels, 0..1, 0..1, 0..1],0));
        } else {
            this.multipoleTranslations = null;
        }

        precomputeTranslations();
    }
    
    public def calculateEnergy() : Double {
        multipoleLowestLevel();
        combineMultipoles();
        transformToLocal();

        return getEnergy();
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        Console.OUT.println("multipole lowest level");
        finish {
            for ((i) in 0..atoms.length-1) {
                val atom = atoms(i);
                val boxLocation = getLowestLevelBoxLocation(atom);
                val boxIndex = FmmBox.getBoxIndex(boxLocation, numLevels);
                async (boxes.dist(boxIndex, numLevels)) {
                    val parentBox = getParentBox(boxIndex, numLevels);
                    val remoteAtom = new MMAtom(atom);
                    var box : FmmBox;
                    atomic {
                        box = boxes(boxIndex, numLevels);
                        if (box == null) {
                            box = new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox);
                            boxes(boxIndex, numLevels) = box;
                        }
                        val leafBox = box as FmmLeafBox!;
                        leafBox.atoms.add(remoteAtom);
                        val boxCentre = leafBox.getCentre(size).sub(remoteAtom.centre);
                        leafBox.multipoleExp.add(MultipoleExpansion.getOlm(remoteAtom.charge, boxCentre, numTerms));
                    }
                }
            }
        }
    }

    /** 
     * Starting at the bottom level, combines multipole expansions for <= 8 child
     * boxes into a single multipole expansion for the parent box.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelRegion : Region(2) = [0..((Math.pow(8,level) as Int)-1),level..level];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            // TODO XTENLANG-1143
            finish ateach ((p1) in Dist.makeUnique(thisLevelDist.places())) {
                finish foreach ((boxIndex1,level) in thisLevelDist | here) {
                    if (boxes(boxIndex1, level) != null) {
                        val child = boxes(boxIndex1, level) as FmmBox!;
                        val childLoc = child.gridLoc;
                        val childExp = child.multipoleExp;
                        val parent = child.parent;
                        at (parent) {
                            val shift = multipoleTranslations(Point.make([here.id, level, (childLoc.x+1)%2, (childLoc.y+1)%2, (childLoc.z+1)%2])) as MultipoleExpansion!;
                            parent.multipoleExp.translateAndAddMultipole(shift, childExp);
                        }
                    }
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
        Console.OUT.println("transform level 2");
        val level2Region : Region(2) = [0..63,2..2];
        val level2Dist = boxes.dist | level2Region;
        // TODO XTENLANG-1143
        finish ateach ((p1) in Dist.makeUnique(level2Dist.places())) {
            finish foreach ((boxIndex1,level) in level2Dist | here) {
                val box1 = boxes(boxIndex1,level) as FmmBox!;
                if (box1 != null) {
                    for ((boxIndex2,level) in level2Dist) {
                        if (boxIndex2 != boxIndex1) {
                            val box2Loc = FmmBox.getBoxLocation(boxIndex2,level);
                            if (box1.wellSeparated(ws, box2Loc)) {
                                val box2MultipoleExp = getMultipoleExpansionLocalCopy(boxIndex2,level);
                                if (box2MultipoleExp != null) {
                                    val translation = box1.getTranslationIndex(box2Loc);
                                    val transform21 = multipoleTransforms(Point.make([here.id, level, -translation.x, -translation.y, -translation.z])) as LocalExpansion!;
                                    box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        for ((thisLevel) in 3..numLevels) {
            Console.OUT.println("transform level " + thisLevel);
            val thisLevelRegion : Region(2) = [0..(Math.pow(8,thisLevel) as Int)-1,thisLevel..thisLevel];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            // TODO XTENLANG-1143
            finish ateach ((p1) in Dist.makeUnique(thisLevelDist.places())) {
                finish foreach ((boxIndex1,level) in thisLevelDist | here) {
                    val box1 = boxes(boxIndex1,level) as FmmBox!;
                    if (box1 != null) {
                        for ((boxIndex2,level) in thisLevelDist) {
                            if (boxIndex2 != boxIndex1) {
                                val box2Loc = FmmBox.getBoxLocation(boxIndex2,level);
                                if (box1.wellSeparated(ws, box2Loc)) {
                                    val parentIndex = getParentIndex(boxIndex2,level);
                                    val box2ParentLoc = FmmBox.getBoxLocation(parentIndex,level-1);
                                    if (!box1.parent.wellSeparated(ws, box2ParentLoc)) {
                                        val box2MultipoleExp = getMultipoleExpansionLocalCopy(boxIndex2,level);
                                        if (box2MultipoleExp != null) {
                                            val translation = box1.getTranslationIndex(box2Loc);
                                            val translateP = Point.make([here.id, level, -translation.x, -translation.y, -translation.z]);
                                            val transform21 = multipoleTransforms(translateP) as LocalExpansion!;
                                            box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                        }
                                    }
                                }
                            }
                        }
                        val box1ParentExp = box1.parent.getLocalExpansionLocalCopy(numTerms);
                        val shift = multipoleTranslations(Point.make([here.id, level, box1.gridLoc.x%2, box1.gridLoc.y%2, box1.gridLoc.z%2])) as MultipoleExpansion!;
                        box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                    }
                }
            }
        }
    }

    def getEnergy() : Double {
        Console.OUT.println("getEnergy");
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
        val lowestLevelRegion : Region(2) = [0..(Math.pow(8,numLevels) as Int)-1,numLevels..numLevels];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        // TODO XTENLANG-1143
        finish ateach ((p1) in Dist.makeUnique(lowestLevelDist.places())) {
            finish foreach ((boxIndex1,level) in lowestLevelDist | here) {
                val box1 = boxes(boxIndex1, numLevels) as FmmLeafBox!;
                if (box1 != null) {
                    // TODO use shared var - XTENLANG-404
                    var thisBoxEnergy : Double = 0.0;
                    val length = box1.atoms.length();
                    for ((atomIndex1) in 0..length-1) {
                        val atom1 = box1.atoms(atomIndex1);
                        val box1Centre = atom1.centre.sub(box1.getCentre(size));
                        val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                        thisBoxEnergy += farFieldEnergy;

                        // direct calculation with all atoms in same box
                        for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                            val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                            val pairEnergy : Double = atom1.pairEnergy(sameBoxAtom);
                            thisBoxEnergy += 2 * pairEnergy;
                        }
                    }
                    if (boxIndex1 > 0) { // EmptyRegion problem
                        // direct calculation with all atoms in non-well-separated boxes
                        val otherBoxRegion : Region(2) = [0..boxIndex1-1,numLevels..numLevels];
                        val otherBoxDist = boxes.dist | otherBoxRegion;
                        for ((boxIndex2,level) in otherBoxDist) {
                            var directBox2Energy : Double = 0.0;
                            val box2Loc = FmmBox.getBoxLocation(boxIndex2,level);
                            if (!box1.wellSeparated(ws, box2Loc)) {
                                val packedAtoms = at(boxes.dist(boxIndex2, level)) {getPackedAtomsForBox(boxIndex2, level)};
                                if (packedAtoms != null) {
                                    for (var i : Int = 0; i < packedAtoms.length(); i+=4) {
                                        val atom2Centre = new Point3d(packedAtoms(i+1), packedAtoms(i+2), packedAtoms(i+3));
                                        val atom2Charge = packedAtoms(i);
                                        for ((atomIndex1) in 0..length-1) {
                                            val atom1 = box1.atoms(atomIndex1);
                                            directBox2Energy += atom1.charge * atom2Charge / atom1.centre.distance(atom2Centre);
                                        }
                                    }
                                    thisBoxEnergy += 2 * directBox2Energy;
                                }
                            }
                        }
                    }
                    val thisBoxEnergyFinal = thisBoxEnergy;
                    async (this) {atomic {fmmEnergy += thisBoxEnergyFinal;}}
                }
            }
        }

        return fmmEnergy;
    }

    /**
     * Precomputes multipole translations and multipole-to-local transformations,
     * which are duplicated at each place.
     */
    private def precomputeTranslations() {
        finish {
            if (numLevels >= 3) {
                ateach ((p) : Point in Dist.makeUnique(Place.places)) {
                    for (val(placeId,level,i,j,k) in multipoleTranslations.dist | here) {
                        dim : Int = Math.pow2(level);
                        sideLength : Double = size / dim;
                        val translationVector = new Vector3d((i*2-1) * 0.5 * sideLength,
                                                         (j*2-1) * 0.5 * sideLength,
                                                         (k*2-1) * 0.5 * sideLength);
                        multipoleTranslations(Point.make([placeId, level, i, j, k])) = MultipoleExpansion.getOlm(translationVector as Tuple3d, numTerms);
                    } 
                }
            }

            ateach ((p) : Point in Dist.makeUnique(Place.places)) {
                for (val(placeId,level,i,j,k) in multipoleTransforms.dist | here) {
                    dim : Int = Math.pow2(level);
                    sideLength : Double = size / dim;
                    val translationVector = new Vector3d(i * sideLength,
                                                     j * sideLength,
                                                     k * sideLength);
                    multipoleTransforms(Point.make([placeId, level, i, j, k])) = LocalExpansion.getMlm(translationVector, numTerms);
                }
            }
        }
    }

    private global def getLowestLevelBoxLocation(atom : MMAtom!) : GridLocation {
        return  GridLocation(atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int);
    }

    private global def getParentBox(childIndex : Int, childLevel : Int) : FmmBox {
        if (childLevel == 2)
            return null;
        val parentIndex = getParentIndex(childIndex, childLevel);
        val parentLevel = childLevel - 1;
        var parent : FmmBox = at (boxes.dist(parentIndex, parentLevel))  boxes(parentIndex, parentLevel);
        if (parent == null) {
            parent = at (boxes.dist(parentIndex, parentLevel)) {createNewParent(parentIndex, parentLevel)};
        }
        return parent;
    }

    private global def createNewParent(parentIndex : int, parentLevel : Int) : FmmBox {
        val grandparent = getParentBox(parentIndex, parentLevel);
        val parentLocation = FmmBox.getBoxLocation(parentIndex, parentLevel);
        val newParent = new FmmBox(parentLevel, parentLocation, numTerms, grandparent);
        boxes(parentIndex, parentLevel) = newParent;
        return newParent;
    }

    private global def getParentIndex(childIndex : Int, childLevel : Int) : Int {
        parentDim : Int = Math.pow2(childLevel-1) as Int;
        val gridLoc = FmmBox.getBoxLocation(childIndex, childLevel);
        return gridLoc.x / 2 * parentDim * parentDim + gridLoc.y / 2 * parentDim + gridLoc.z / 2;
    }

    /**
     * Gets the pairwise interaction energy between the given atom and all atoms in the given box.
     */
    private global def getPairwiseInteractionForBox(atom : MMAtom, boxIndex : Int) : Double {
        var boxEnergy : Double = 0.0;
        val box = boxes(boxIndex, numLevels) as FmmLeafBox!;
        if (box != null) {
            for ((atomIndex2) in 0..box.atoms.length()-1) {
                val atom2 = box.atoms(atomIndex2) as MMAtom!;
                val pairEnergy : Double = atom2.pairEnergy(atom);
                boxEnergy += 2 * pairEnergy;
            }
        }
        return boxEnergy;
    }

    private global def getMultipoleForBox(boxIndex : Int, level : Int) : MultipoleExpansion {
        val box = boxes(boxIndex, level) as FmmBox!;
        return box == null ? null : box.multipoleExp;
    }

    private global def getPackedAtomsForBox(boxIndex : Int, level : Int) {
        val box = boxes(boxIndex,level) as FmmLeafBox!;
        if (box != null) {
            return box.getPackedAtoms();
        } else {
            return null;
        }
    }

    /**
     * TODO this is a workaround due to lack of Array copy facility - XTENLANG-787
     * @return a local copy at the current place of this box's multipole expansion
     */
    public global def getMultipoleExpansionLocalCopy(boxIndex : Int, level : Int) : MultipoleExpansion! {
        val data = at (boxes.dist(boxIndex, level)) {boxes(boxIndex, level) != null? Expansion.getData(numTerms, (boxes(boxIndex, level) as FmmBox!).multipoleExp) : null};
        if (data != null) {
            return new MultipoleExpansion(numTerms, data);
        } else {
            return null;
        }
    }
}

