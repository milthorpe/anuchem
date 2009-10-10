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
    val boxes : Array[Box[FmmBox]]{rank==2};

    val atoms : ValRail[Atom];

    /** A cache of transformations from multipole to local at the same level. */
    val multipoleTransforms : Array[LocalExpansion]{rank==5};

    /** A cache of multipole translations between parent box centres and child box centres. */
    val multipoleTranslations : Array[MultipoleExpansion]{rank==5};

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


        var boxRegion : Region{rank==2} = [0..63, 2..2];
        for ((i) in 3..numLevels) {
            rNextLevel : Region{rank==2} = [0..(Math.pow(8,i) as Int)-1, i..i];
            boxRegion = boxRegion || rNextLevel;
        }
        Console.OUT.println("boxes: " + boxRegion);

        val boxDistribution : Dist{rank==2} = Dist.makeBlock(boxRegion, 0);
        Console.OUT.println("dist: " + boxDistribution);
        
        // all boxes are null to start.  they will be initialised as needed.
        this.boxes = Array.make[Box[FmmBox]](boxDistribution);

        // two special arrays distributed to all places (this is done by replicating the first index in a cyclic dist across Place.PLACES)
        var wellSpacedLimit : Region(5) = [0..Place.MAX_PLACES-1,2..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransformRegion : Region(5) = wellSpacedLimit - ([0..Place.MAX_PLACES-1,2..numLevels,-ws..ws,-ws..ws,-ws..ws] as Region);
        this.multipoleTransforms = Array.make[LocalExpansion](Dist.makeCyclic(multipoleTransformRegion,0));
        this.multipoleTranslations = Array.make[MultipoleExpansion](Dist.makeCyclic([0..Place.MAX_PLACES-1,2..numLevels, 0..1, 0..1, 0..1],0));
    }
    
    public def calculateEnergy() : Double {
        // TODO XTENLANG-513
        val numTerms = this.numTerms;
        val size = this.size;
        val multipoleTranslations = this.multipoleTranslations;

        // precompute multipole translations
        val d: Dist = Dist.makeUnique(Place.places);
        finish ateach ((p) : Point in d) {
            val myPortion = multipoleTranslations.dist | here;
            foreach (val(placeId,level,i,j,k) in multipoleTranslations.dist | here) {
                dim : Int = Math.pow2(level);
                sideLength : Double = size / dim;
                translationVector : Vector3d = new Vector3d((i*2-1) * 0.5 * sideLength,
                                                                 (j*2-1) * 0.5 * sideLength,
                                                                 (k*2-1) * 0.5 * sideLength);
                multipoleTranslations(Point.make([placeId, level, i, j, k])) = MultipoleExpansion.getOlm(translationVector, numTerms);
            } 
        }

        Console.OUT.println("multipoleLowestLevel");

        multipoleLowestLevel();
        /*
        var nonEmpty : Int = 0;
        for (val (i,j) in boxes.region) {
            if (boxes(i,j) != null) {
                nonEmpty++;
                Console.OUT.println("boxes(" + i + "," + j + ") multipole = " + boxes(i,j).value.multipoleExp);
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
        
        for ((i) in 0..atoms.length-1) {
            val atom = atoms(i);
            val boxLocation = getLowestLevelBoxLocation(atom);
            val boxIndex = FmmBox.getBoxIndex(boxLocation, numLevels);
            //Console.OUT.println("atoms(" + i + ") => box(" + boxIndex + ")");
            val parentBox = getParentBox(boxIndex, numLevels);
            //Console.OUT.println("boxIndex = " + boxIndex + " numLevels = " + numLevels);

            // TODO XTENLANG-513
            val numLevels = this.numLevels;
            val numTerms = this.numTerms;
            val boxes = this.boxes;
            val size = this.size;
            at (boxes.dist(boxIndex, numLevels)) {
                var box : Box[FmmBox] = boxes(boxIndex, numLevels);
                if (box == null) {
                    box = new Box[FmmBox](new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox));
                    boxes(boxIndex, numLevels) = box;
                }
                val remoteAtom = new Atom(atom);
                val leafBox = box.value as FmmLeafBox!;
                leafBox.atoms.add(remoteAtom);
                val boxCentre = leafBox.getCentre(size).sub(remoteAtom.centre);
                leafBox.multipoleExp.add(MultipoleExpansion.getOlm(remoteAtom.charge, boxCentre, numTerms));
            }
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            val thisLevelRegion :Region{rank==2} = [0..((Math.pow(8,level) as Int)-1),level..level];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach ((childIndex,childLevel) in thisLevelDist) {
                //Console.OUT.println("childIndex = " + childIndex + " childLevel = " + childLevel + " dist " + boxes.dist(childIndex,childLevel));
                val bChild = boxes(childIndex,childLevel);
                if (bChild != null) {
                    val child = bChild.value as FmmBox!;
                    val bParent = child.parent;
                    val parentExp = at(bParent.value) {bParent.value.multipoleExp};
                    val shift = multipoleTranslations(Point.make([here.id, childLevel, (child.gridLoc(0)+1)%2, (child.gridLoc(1)+1)%2, (child.gridLoc(2)+1)%2]));
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
        // TODO XTENLANG-513
        val ws = this.ws;
        val numLevels = this.numLevels;
        val boxes = this.boxes;
        val numTerms = this.numTerms;
        val size = this.size;
        val multipoleTransforms = this.multipoleTransforms;
        val multipoleTranslations = this.multipoleTranslations;

        // precompute child translations
        val d: Dist = Dist.makeUnique(Place.places);
        finish ateach ((p) : Point in d) {
            foreach (val(placeId, level,i,j,k) in multipoleTransforms.dist | here) {
                dim : Int = Math.pow2(level);
                sideLength : Double = size / dim;
                translationVector : Vector3d = new Vector3d(i * sideLength,
                                                            j * sideLength,
                                                            k * sideLength);
                multipoleTransforms(Point.make([placeId, level, i, j, k])) = LocalExpansion.getMlm(translationVector, numTerms);
            }
        } 

        Console.OUT.println("level 2");

        val level2Region : Region{rank==2} = [0..63,2..2];
        val level2Dist = boxes.dist | level2Region;
        finish ateach (p1(boxIndex1,level) in level2Dist) {
            val bBox1 = boxes(p1);
            if (bBox1 != null) {
                val box1 = bBox1.value as FmmBox!;
                val box1MultipoleExp = box1.multipoleExp;
                //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + 2 + ")");
                for ((boxIndex2) in 0..boxIndex1-1) { 
                    async (boxes.dist(boxIndex2,level)) {
                        val bBox2 = boxes(boxIndex2, level);
                        if (bBox2 != null) {
                            val box2 = bBox2.value as FmmBox!;
                            //Console.OUT.println("... and box(" + 2 + "," + boxIndex2 + ")");
                            if (box2.wellSeparated(ws, box1)) {
                                val translation = box2.getTranslationIndex(box1);
                                val transform12 = multipoleTransforms(Point.make([here.id, level, -translation(0), -translation(1), -translation(2)]));
                                box2.localExp.transformAndAddToLocalDist(transform12, box1MultipoleExp);
                                val box2MultipoleExp = box2.multipoleExp;
                                at (bBox1.location) {
                                    val transform21 = multipoleTransforms(Point.make([here.id, level, translation(0), translation(1), translation(2)]));
                                    box1.localExp.transformAndAddToLocalDist(transform21, box2MultipoleExp);
                                }
                            }
                        }
                    }
                }
            }
        }

        Console.OUT.println("remaining levels");
        
        for ((thisLevel) in 3..numLevels) {
            val thisLevelRegion : Region{rank==2} = [0..(Math.pow(8,thisLevel) as Int)-1,thisLevel..thisLevel];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach ((boxIndex1,level) in thisLevelDist) {
                val bBox1 = boxes(boxIndex1, level);
                if (bBox1 != null) {
                    val box1 = bBox1.value as FmmBox!;
                    val box1MultipoleExp = box1.multipoleExp;
                    val box1Parent = at(box1.parent.location) {box1.parent.value};
                    //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ")");
                    for ((boxIndex2) in 0..boxIndex1-1) {
                        async (boxes.dist(boxIndex2,level)) {
                            val bBox2 = boxes(boxIndex2, level);
                            if (bBox2 != null) {
                                val box2 = bBox2.value as FmmBox!;
                                val box2Parent = at(box2.parent.location) {box2.parent.value};
                                //Console.OUT.println("... and box(" + level + "," + boxIndex2 + ")");
                                if (!(at (box2.parent.location) {box2.parent.value.wellSeparated(ws, box1Parent)})) {
                                    if (box2.wellSeparated(ws, box1)) {
                                        val translation = box2.getTranslationIndex(box1);
                                        val transform12 = multipoleTransforms(Point.make([here.id, level, -translation(0), -translation(1), -translation(2)]));
                                        box2.localExp.transformAndAddToLocalDist(transform12, box1MultipoleExp);
                                        val box2MultipoleExp = box2.multipoleExp;
                                        at (bBox1.location) {
                                            val transform21 = multipoleTransforms(Point.make([here.id, level, translation(0), translation(1), translation(2)]));
                                            box1.localExp.transformAndAddToLocalDist(transform21, box2MultipoleExp);
                                        }   
                                    }
                                }
                            }
                        }
                    }
                    val shift = multipoleTranslations(Point.make([here.id, level, box1.gridLoc(0)%2, box1.gridLoc(1)%2, box1.gridLoc(2)%2]));
                    box1.localExp.translateAndAddLocal(shift, box1Parent.localExp);
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
                val pairEnergy : Double = pairEnergy(atoms(i), atoms(j));
                directEnergy += 2 * pairEnergy;
            }
        }
        Console.OUT.println("directEnergy = " + directEnergy);
        */

        // TODO XTENLANG-513
        val ws = this.ws;
        val numLevels = this.numLevels;
        val size = this.size;
        val boxes = this.boxes;

        val lowestLevelDist = boxes.dist | [0..(Math.pow(8,numLevels) as Int)-1,numLevels..numLevels];
        finish ateach ((boxIndex1,level) in lowestLevelDist) {
            //Console.OUT.println("boxIndex1 = " + boxIndex1);
            val bBox1 = boxes(boxIndex1, numLevels);
            if (bBox1 != null) {
                val box1 = bBox1.value as FmmLeafBox!;
                //Console.OUT.println("getEnergy: box(" + boxIndex1 + "," + numLevels + ")");
                for ((atomIndex1) in 0..box1.atoms.length()-1) {
                    // TODO should be able to use a shared var for atom energy
                    val atom1 = box1.atoms(atomIndex1);
                    val box1Centre = atom1.centre.sub(box1.getCentre(size));
                    val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                    //Console.OUT.println("farFieldEnergy " + farFieldEnergy + " at " + this.location);
                    at (this.location){atomic {fmmEnergy += farFieldEnergy;}}

                    //Console.OUT.println("direct - same box");
                    // direct calculation with all atoms in same box
                    for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                        val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                        val pairEnergy : Double = atom1.pairEnergy(sameBoxAtom);
                        at (this.location){atomic {fmmEnergy += 2 * pairEnergy;}}
                    }

                    //Console.OUT.println("direct - non-well-sep");
                    // direct calculation with all atoms in non-well-separated boxes
                    val otherBoxDist = boxes.dist | [0..boxIndex1-1,numLevels..numLevels];
                    finish ateach ((boxIndex2,level) in otherBoxDist) {
                        val bBox2 = boxes(boxIndex2, level);
                        if (bBox2 != null) {
                            val box2 = bBox2.value as FmmLeafBox;
                            //Console.OUT.println(boxIndex1 + " vs. " + boxIndex2 + " ws = " + ws);
                            if (!box2.wellSeparated(ws, box1)) {
                                for ((atomIndex2) in 0..box2.atoms.length()-1) {
                                    //Console.OUT.println("pair energy: " + boxIndex1 + "-" + atomIndex1 + " " + boxIndex2 + "-" + atomIndex2);
                                    val atom2 = box2.atoms(atomIndex2);
                                    val pairEnergy : Double = atom2.pairEnergy(atom1);
                                    //Console.OUT.println("pairEnergy " + here);
                                    at (this.location){atomic {fmmEnergy += 2 * pairEnergy;}}
                                }
                            }
                        }
                    }
                }
            }
        }

        return fmmEnergy;
    }

    def getLowestLevelBoxLocation(atom : Atom) : ValRail[Int]{length==3} {
        index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        return index;
    }

    def getParentBox(childIndex : Int, childLevel : Int) : Box[FmmBox] {
        //Console.OUT.println("getParentBox(" + childIndex + ", " + childLevel + ")");
        if (childLevel == 2)
            return null;
        parentIndex : Int = getParentIndex(childIndex, childLevel);
        parentLevel : Int = childLevel - 1;
        var parent : Box[FmmBox] = at (boxes.dist(parentIndex, parentLevel)) {boxes(parentIndex, parentLevel)};
        if (parent == null) {
            parent = at (boxes.dist(parentIndex, parentLevel)) {getNewParent(parentIndex, parentLevel)};
        }
        return parent;
    }

    def getNewParent(parentIndex : int, parentLevel : Int) : FmmBox {
        val grandparent = getParentBox(parentIndex, parentLevel);
        val parentLocation = getBoxLocation(parentIndex, parentLevel);
        val newParent = new FmmBox(parentLevel, parentLocation, numTerms, grandparent);
        if (grandparent != null)
            boxes(parentIndex, parentLevel) = newParent;
        return newParent;
    }

    def getBoxLocation(index : Int, level : Int) : ValRail[Int]{length==3} {
        dim : Int = Math.pow2(level) as Int;
        gridLoc : ValRail[Int]{length==3} = [ index / (dim * dim), (index / dim) % dim, index % dim ];
        return gridLoc;
    }

    def getParentIndex(childIndex : Int, childLevel : Int) : Int {
        parentDim : Int = Math.pow2(childLevel-1) as Int;
        gridLoc : ValRail[Int] = getBoxLocation(childIndex, childLevel);
        return gridLoc(0) / 2 * parentDim * parentDim + gridLoc(1) / 2 * parentDim + gridLoc(2) / 2;
    }
}

