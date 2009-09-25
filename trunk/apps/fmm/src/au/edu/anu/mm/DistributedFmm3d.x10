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

    /** The maximum number of boxes in the tree, if there are no empty boxes at the lowest level. */
    val maxBoxes : Int;

    /** All boxes in the octree division of space. */
    val boxes : Array[FmmBox]{rank==2};

    val atoms : ValRail[Atom];

    /** A cache of transformations from multipole to local at the same level. */
    val multipoleTransforms : Array[LocalExpansion]{rank==4};

    /** A cache of multipole translations between parent box centres and child box centres. */
    val multipoleTranslations : Array[MultipoleExpansion]{rank==4};

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

        var nBox : Int = 1;
        for ((i) in 2..numLevels) {
            nBox += Math.pow(8,i) as Int;
        }
        this.maxBoxes = nBox;
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
        this.boxes = Array.make[FmmBox](boxDistribution);

    
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
                var box : FmmLeafBox = boxes(boxIndex, numLevels) as FmmLeafBox;
                if (box == null) {
                    val newBox = new FmmLeafBox(numLevels, boxLocation, numTerms, parentBox);
                    boxes(boxIndex, numLevels) = newBox;
                    box = newBox;
                }
                val remoteAtom = new Atom(atom);
                box.addAtom(remoteAtom);
                v : Tuple3d = box.getCentre(size).sub(remoteAtom.centre);
                olm : MultipoleExpansion = MultipoleExpansion.getOlm(remoteAtom.charge, v, numTerms);
                box.multipoleExp.add(olm);
            }
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def combineMultipoles() {
        for (var level: Int = numLevels; level > 2; level--) {
            val thisLevelDist = boxes.dist | [0..((Math.pow(8,level) as Int)-1),level..level];
            finish ateach ((childIndex,childLevel) in thisLevelDist) {
            //for ((childIndex) in 0..(Math.pow(8,level) as Int)-1) {
            //    childLevel : Int = level;
                //Console.OUT.println("childIndex = " + childIndex + " childLevel = " + childLevel + " dist " + boxes.dist(childIndex,childLevel));
            //    at (boxes.dist(childIndex,childLevel)) {
                    val child = boxes(childIndex,childLevel);
                    if (child != null) {
                        val parent = child.parent;
                        val shift = at (Place.FIRST_PLACE) {multipoleTranslations(childLevel, (child.gridLoc(0)+1)%2, (child.gridLoc(1)+1)%2, (child.gridLoc(2)+1)%2)};
                        MultipoleExpansion.translateAndAddMultipole(shift, child.multipoleExp, parent.multipoleExp);
                    }
                }
            //}
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
        finish foreach (val(level,i,j,k) in multipoleTransforms.region) {
            dim : Int = Math.pow2(level);
            sideLength : Double = size / dim;
            translationVector : Vector3d = new Vector3d(i * sideLength,
                                                        j * sideLength,
                                                        k * sideLength);
            multipoleTransforms(level, i, j, k) = LocalExpansion.getMlm(translationVector, numTerms);
        } 

        Console.OUT.println("level 2");

        // TODO XTENLANG-513
        val ws = this.ws;
        val numLevels = this.numLevels;
        val boxes = this.boxes;

        val level2Dist = boxes.dist | [0..63,2..2];
        finish ateach ((boxIndex1,level) in level2Dist) {
            val box1 = boxes(boxIndex1, level);
            if (box1 != null) {
                //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + 2 + ")");
                for ((boxIndex2) in 0..boxIndex1-1) { 
                    async (boxes.dist(boxIndex2,level)) {
                        val box2 = boxes(boxIndex2, level);
                        if (box2 != null) {
                            //Console.OUT.println("... and box(" + 2 + "," + boxIndex2 + ")");
                            if (box2.wellSeparated(ws, box1)) {
                                val translation = box2.getTranslationIndex(box1);
                                val transform12 = at (Place.FIRST_PLACE) {multipoleTransforms(level, -translation(0), -translation(1), -translation(2))};
                                box2.localExp.transformAndAddToLocal(transform12, at (boxes.dist(boxIndex1,2)) {box1.multipoleExp});
                                val box2MultipoleExp = at (boxes.dist(boxIndex2,2)) {box2.multipoleExp};
                                at (boxes.dist(boxIndex1,level)) {
                                    val transform21 = at (Place.FIRST_PLACE) {multipoleTransforms(level, translation(0), translation(1), translation(2))};
                                    box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                }
                            }
                        }
                    }
                }
            }
        }

        Console.OUT.println("remaining levels");

        
        for ((thisLevel) in 3..numLevels) {
            val thisLevelDist = boxes.dist | [0..(Math.pow(8,thisLevel) as Int)-1,thisLevel..thisLevel];
            finish ateach ((boxIndex1,level) in thisLevelDist) {
                    val box1 = boxes(boxIndex1, level);
                    if (box1 != null) {
                        //Console.OUT.println("transformToLocal: box(" + boxIndex1 + "," + level + ")");
                        for ((boxIndex2) in 0..boxIndex1-1) {
                            async (boxes.dist(boxIndex2,level)) {
                                box2 : FmmBox = boxes(boxIndex2, level);
                                if (box2 != null) {
                                    //Console.OUT.println("... and box(" + level + "," + boxIndex2 + ")");
                                    if (!(at (box2.parent.location) {box2.parent.wellSeparated(ws, box1.parent)})) {
                                        if (box2.wellSeparated(ws, box1)) {
                                            val translation = box2.getTranslationIndex(box1);
                                            transform12 : LocalExpansion = at (Place.FIRST_PLACE){multipoleTransforms(level, -translation(0), -translation(1), -translation(2))};
                                            box2.localExp.transformAndAddToLocal(transform12, at (boxes.dist(boxIndex1,level)) {box1.multipoleExp});
                                            val box2MultipoleExp = at (boxes.dist(boxIndex2,2)) {box2.multipoleExp};
                                            at (boxes.dist(boxIndex1,2)) {
                                                transform21 : LocalExpansion = at (Place.FIRST_PLACE){multipoleTransforms(level, translation(0), translation(1), translation(2))};
                                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                            }   
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    shift : MultipoleExpansion = multipoleTranslations(level, box1.gridLoc(0)%2, box1.gridLoc(1)%2, box1.gridLoc(2)%2);
                    LocalExpansion.translateAndAddLocal(shift, box1.parent.localExp, box1.localExp);
                }
            //}
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
        //for ((boxIndex1) in 0..(Math.pow(8,numLevels) as Int)-1) { 
            //Console.OUT.println("boxIndex1 = " + boxIndex1);
        //    at (boxes.dist(boxIndex1,numLevels)) {
                val box1 = boxes(boxIndex1, numLevels) as FmmLeafBox;
                if (box1 != null) {
                    //Console.OUT.println("getEnergy: box(" + boxIndex1 + "," + numLevels + ")");
                    for ((atomIndex1) in 0..box1.atoms.length()-1) {
                        atom1 : Atom = box1.atoms(atomIndex1);
                        v : Tuple3d = atom1.centre.sub(box1.getCentre(size));
                        farFieldEnergy : Double = box1.localExp.getPotential(atom1.charge, v);
                        //Console.OUT.println("farFieldEnergy " + farFieldEnergy + " at " + this.location);
                        at (this.location){atomic {fmmEnergy += farFieldEnergy;}}

                        //Console.OUT.println("direct - same box");
                        // direct calculation with all atoms in same box
                        for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                            sameBoxAtom : Atom = box1.atoms(sameBoxAtomIndex);
                            val pairEnergy : Double = atom1.pairEnergy(sameBoxAtom);
                            at (this.location){atomic {fmmEnergy += 2 * pairEnergy;}}
                        }

                        //Console.OUT.println("direct - non-well-sep");
                        // direct calculation with all atoms in non-well-separated boxes
                        finish foreach ((boxIndex2) in 0..boxIndex1-1) {
                            at (boxes.dist(boxIndex2,numLevels)) {
                                box2 : FmmLeafBox = boxes(boxIndex2, numLevels) as FmmLeafBox;
                                if (box2 != null) {
                                    //Console.OUT.println(boxIndex1 + " vs. " + boxIndex2 + " ws = " + ws);
                                    val box1 = at (boxes.dist(boxIndex1,numLevels)) {boxes(boxIndex1, numLevels)};
                                    if (!(box2.wellSeparated(ws, box1))) {
                                        for ((atomIndex2) in 0..box2.atoms.length()-1) {
                                            //Console.OUT.println("pair energy: " + boxIndex1 + "-" + atomIndex1 + " " + boxIndex2 + "-" + atomIndex2);
                                            atom2 : Atom = box2.atoms(atomIndex2);
                                            val pairEnergy : Double = atom2.pairEnergy(atom1);
                                            //Console.OUT.println("pairEnergy " + here);
                                            at (this.location) {atomic {fmmEnergy += 2 * pairEnergy;}}
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            //}
        }

        return fmmEnergy;
    }

    def getLowestLevelBoxLocation(atom : Atom) : ValRail[Int]{length==3} {
        index : ValRail[Int]{length==3} = [ atom.centre.i / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.j / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int, atom.centre.k / size * dimLowestLevelBoxes + dimLowestLevelBoxes / 2 as Int ];
        return index;
    }

    def getParentBox(childIndex : Int, childLevel : Int) : FmmParentBox {
        //Console.OUT.println("getParentBox(" + childIndex + ", " + childLevel + ")");
        if (childLevel == 2)
            return null;
        parentIndex : Int = getParentIndex(childIndex, childLevel);
        parentLevel : Int = childLevel - 1;
        var parent : FmmParentBox = at (boxes.dist(parentIndex, parentLevel)) {
            boxes(parentIndex, parentLevel) as FmmParentBox
        };
        
        if (parent == null) {
            parent = at (boxes.dist(parentIndex, parentLevel)) {getNewParent(parentIndex, parentLevel)};
        }
        return parent;
    }

    def getNewParent(parentIndex : int, parentLevel : Int) : FmmParentBox {
        grandparent : FmmParentBox = getParentBox(parentIndex, parentLevel);
        parentLocation : ValRail[Int]{length==3} = getBoxLocation(parentIndex, parentLevel);
        newParent : FmmParentBox = new FmmParentBox(parentLevel, parentLocation, numTerms, grandparent);
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

