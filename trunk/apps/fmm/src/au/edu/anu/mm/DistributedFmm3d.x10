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
public class DistributedFmm3d {
    /** The number of levels in the octree. */
    public global val numLevels : Int;

    /** The number of lowest level boxes along one side of the (cubic) 3D space. */
    public global val dimLowestLevelBoxes : Int;

    /** The cartesian location of the top-left-front corner of the simulation cube. */
    public val topLeftFront : Point3d;

    /** The length of a side of the simulation cube. */
    public global val size : Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public global val numTerms : Int;

    /** The well-separatedness parameter ws. */
    public global val ws : Int;

    /** All boxes in the octree division of space. */
    global val boxes : Array[FmmBox](2);

    val atoms : ValRail[Atom];

    /** A cache of transformations from multipole to local at the same level. */
    global val multipoleTransforms : Array[LocalExpansion](5);

    /** A cache of multipole translations between parent box centres and child box centres. */
    global val multipoleTranslations : Array[MultipoleExpansion](5);

    // TODO use shared local variable within getEnergy() - 
    // not currently possible due to <a href="http://jira.codehaus.org/browse/XTENLANG-505"/>
    var fmmEnergy : Double = 0.0;

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

        var boxRegion : Region(2) = [0..63, 2..2];
        var boxDistribution : Dist(2) = Dist.makeBlock(boxRegion, 0);
        for ((i) in 3..numLevels) {
            val nextLevelRegion : Region(2) = [0..(Math.pow(8,i) as Int)-1, i..i];
            val nextLevelDist = Dist.makeBlock(nextLevelRegion, 0);
            boxDistribution = boxDistribution || nextLevelDist;
        }
        Console.OUT.println("dist: " + boxDistribution);
        
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[FmmBox](boxDistribution);

        // two special arrays distributed to all places (this is done by replicating the first index in a cyclic dist across Place.PLACES)
        var wellSpacedLimit : Region(5) = [0..Place.MAX_PLACES-1,2..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransformRegion : Region(5) = wellSpacedLimit - ([0..Place.MAX_PLACES-1,2..numLevels,-ws..ws,-ws..ws,-ws..ws] as Region);
        this.multipoleTransforms = Array.make[LocalExpansion](Dist.makeCyclic(multipoleTransformRegion,0));
        if (numLevels >= 3) {
            this.multipoleTranslations = Array.make[MultipoleExpansion](Dist.makeCyclic([0..Place.MAX_PLACES-1,3..numLevels, 0..1, 0..1, 0..1],0));
        } else {
            this.multipoleTranslations = null;
        }
    }
    
    public def calculateEnergy() : Double {
        if (numLevels >= 3) {
            // precompute multipole translations
            val d: Dist = Dist.makeUnique(Place.places);
            finish ateach ((p) : Point in d) {
                for (val(placeId,level,i,j,k) in multipoleTranslations.dist | here) {
                    dim : Int = Math.pow2(level);
                    sideLength : Double = size / dim;
                    translationVector : Vector3d = new Vector3d((i*2-1) * 0.5 * sideLength,
                                                                     (j*2-1) * 0.5 * sideLength,
                                                                     (k*2-1) * 0.5 * sideLength);
                    multipoleTranslations(Point.make([placeId, level, i, j, k])) = MultipoleExpansion.getOlm(translationVector, numTerms);
                } 
            }
        }

        Console.OUT.println("multipoleLowestLevel");

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

        Console.OUT.println("combine");
        combineMultipoles();
        Console.OUT.println("transform");
        transformToLocal();
        /*
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                Console.OUT.println("boxes(" + i + "," + j + ") local = " + boxes(i,j).localExp);
            }
        }
        */
        Console.OUT.println("getEnergy");
        return getEnergy();
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        finish {
            for ((i) in 0..atoms.length-1) {
                val atom = atoms(i);
                val boxLocation = getLowestLevelBoxLocation(atom);
                val boxIndex = FmmBox.getBoxIndex(boxLocation, numLevels);
                //Console.OUT.println("atoms(" + i + ") => box(" + boxIndex + ")");
                val parentBox = getParentBox(boxIndex, numLevels);
                //Console.OUT.println("boxIndex = " + boxIndex + " numLevels = " + numLevels);

                at (boxes.dist(boxIndex, numLevels)) {
                    atomic {
                        var box : FmmBox = boxes(boxIndex, numLevels);
                        if (box == null) {
                            box = new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox);
                            boxes(boxIndex, numLevels) = box;
                        }
                        val remoteAtom = new Atom(atom);
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
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            val thisLevelRegion : Region(2) = [0..((Math.pow(8,level) as Int)-1),level..level];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            ateach (p in thisLevelDist) {
                //Console.OUT.println("childIndex = " + childIndex + " childLevel = " + childLevel + " dist " + boxes.dist(childIndex,childLevel));
                val child = boxes(p) as FmmBox!;
                if (child != null) {
                    val parent = child.parent as FmmBox!;
                    val shift = multipoleTranslations(Point.make([here.id, child.level, (child.gridLoc.x+1)%2, (child.gridLoc.y+1)%2, (child.gridLoc.z+1)%2])) as MultipoleExpansion!;
                    parent.multipoleExp.translateAndAddMultipole(shift, child.multipoleExp);
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
        val d: Dist = Dist.makeUnique(Place.places);
        finish ateach ((p) : Point in d) {
            for (val(placeId,level,i,j,k) in multipoleTransforms.dist | here) {
                dim : Int = Math.pow2(level);
                sideLength : Double = size / dim;
                translationVector : Vector3d = new Vector3d(i * sideLength,
                                                            j * sideLength,
                                                            k * sideLength);
                multipoleTransforms(Point.make([placeId, level, i, j, k])) = LocalExpansion.getMlm(translationVector, numTerms);
            }
        }

        Console.OUT.println("level 2");

        val level2Region : Region(2) = [0..63,2..2];
        val level2Dist = boxes.dist | level2Region;
        finish ateach (p1(boxIndex1,level) in level2Dist) {
            val box1 = boxes(p1) as FmmBox!;
            if (box1 != null) {
                //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + 2 + ")");
                for (p2(boxIndex2,level) in level2Dist) {
                    if (boxIndex2 != boxIndex1) {
                        val box2Loc = FmmBox.getBoxLocation(boxIndex2,level);
                        if (box1.wellSeparated(ws, box2Loc)) {
                            val box2MultipoleExp = at (boxes.dist(p2)) getMultipoleForBox(p2);
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(box2Loc);
                                val transform21 = multipoleTransforms(Point.make([here.id, level, -translation.x, -translation.y, -translation.z])) as LocalExpansion!;
                                //val box2MultipoleExp = box2MultipoleExpFuture();
                                //Console.OUT.println("added box(" + boxIndex1 + "," + 2 + ") to box(" + boxIndex2 + "," + 2 + ")");
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                    }
                }
            }
        }

        Console.OUT.println("remaining levels");
        
        for ((thisLevel) in 3..numLevels) {
            val thisLevelRegion : Region(2) = [0..(Math.pow(8,thisLevel) as Int)-1,thisLevel..thisLevel];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach (p1(boxIndex1,level) in thisLevelDist) {
                val box1 = boxes(p1) as FmmBox!;
                if (box1 != null) {
                    val box1Parent = box1.parent as FmmBox!;
                    //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ")");
                    for (p2(boxIndex2,level) in thisLevelDist) {
                        if (boxIndex2 != boxIndex1) {
                            val box2Loc = FmmBox.getBoxLocation(boxIndex2,level);
                            if (box1.wellSeparated(ws, box2Loc)) {
                                val parentIndex = getParentIndex(boxIndex2,level);
                                val box2ParentLoc = FmmBox.getBoxLocation(parentIndex,level-1);
                                if (box1.parent.wellSeparated(ws, box2ParentLoc)) {
                                    val box2MultipoleExp = at (boxes.dist(p2)) getMultipoleForBox(p2);
                                    val translation = box1.getTranslationIndex(box2Loc);
                                    val transform21 = multipoleTransforms(Point.make([here.id, level, -translation.x, -translation.y, -translation.z])) as LocalExpansion!;
                                    //Console.OUT.println("added box(" + boxIndex1 + "," + level + ") to box(" + boxIndex2 + "," + level + ")");
                                    box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                }
                            }
                        }
                    }
                    val shift = multipoleTranslations(Point.make([here.id, level, box1.gridLoc.x%2, box1.gridLoc.y%2, box1.gridLoc.z%2])) as MultipoleExpansion!;
                    box1.localExp.translateAndAddLocal(shift, box1Parent.localExp);
                    //Console.OUT.println("added box1(" + boxIndex1 + "," + level + ") parent");
                }
            }
        }
        //Console.OUT.println("wellSep = " + wellSep + " nearField = " + nearField);
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
        val lowestLevelRegion : Region(2) = [0..(Math.pow(8,numLevels) as Int)-1,numLevels..numLevels];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        finish ateach ((boxIndex1,level) in lowestLevelDist) {
            //at (boxes.dist(boxIndex1,level)) {
            //Console.OUT.println("boxIndex1 = " + boxIndex1);
            val box1 = boxes(boxIndex1, numLevels) as FmmLeafBox!;
            if (box1 != null) {
                shared var thisBoxEnergy : Double = 0.0;
                //Console.OUT.println("getEnergy: box(" + boxIndex1 + "," + numLevels + ")");
                for ((atomIndex1) in 0..box1.atoms.length()-1) {
                    // TODO should be able to use a shared var for atom energy
                    val atom1 = box1.atoms(atomIndex1) as Atom!;
                    val box1Centre = atom1.centre.sub(box1.getCentre(size));
                    val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                    //Console.OUT.println("farFieldEnergy " + farFieldEnergy + " at " + this.home);
                    thisBoxEnergy += farFieldEnergy;

                    //Console.OUT.println("direct - same box");
                    // direct calculation with all atoms in same box
                    for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                        val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                        val pairEnergy : Double = atom1.pairEnergy(sameBoxAtom);
                        thisBoxEnergy += 2 * pairEnergy;
                    }

                    if (boxIndex1 > 0) { // EmptyRegion problem
                        //Console.OUT.println("direct - non-well-sep");
                        // direct calculation with all atoms in non-well-separated boxes
                        val otherBoxRegion : Region(2) = [0..boxIndex1-1,numLevels..numLevels];
                        val otherBoxDist = boxes.dist | otherBoxRegion;
                        for ((boxIndex2,level) in otherBoxDist) {
                            val box2 = boxes(boxIndex2, level) as FmmLeafBox!;
                            if (box2 != null) {
                                //Console.OUT.println(boxIndex1 + " vs. " + boxIndex2 + " ws = " + ws);
                                if (!box2.wellSeparated(ws, box1)) {
                                    val box2Energy : Double = at (boxes.dist(boxIndex2, level)) getPairwiseInteractionForBox(atom1, box2 as FmmLeafBox!);
                                    //Console.OUT.println("pair energy: " + boxIndex1 + "-" + atomIndex1 + " with box " + boxIndex2 + " = " + box2Energy);
                                    thisBoxEnergy += box2Energy;
                                }
                            }
                        }
                    }
                }
                val thisBoxEnergyFinal = thisBoxEnergy;
                at (this.home) {atomic {fmmEnergy += thisBoxEnergyFinal;}}
            }
        }

        return fmmEnergy;
    }

    private global def getLowestLevelBoxLocation(atom : Atom) : GridLocation {
        return  GridLocation(atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int);
        /*index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        return index;*/
    }

    private global def getParentBox(childIndex : Int, childLevel : Int) : FmmBox {
        //Console.OUT.println("getParentBox(" + childIndex + ", " + childLevel + ")");
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
    public global def getPairwiseInteractionForBox(atom : Atom, box : FmmLeafBox!) : Double {
        var boxEnergy : Double = 0.0;
        for ((atomIndex2) in 0..box.atoms.length()-1) {
            val atom2 = box.atoms(atomIndex2) as Atom!;
            val pairEnergy : Double = atom2.pairEnergy(atom);
            boxEnergy += 2 * pairEnergy;
        }
        return boxEnergy;
    }

    public global def getMultipoleForBox(p : Point(2)) : MultipoleExpansion {
        val box = boxes(p) as FmmBox!;
        return box == null ? null : box.multipoleExp;
    }
}

