package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

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
    public global val lowestLevelDim : Int;

    /** The side length of the cube. */
    public global val size : Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public global val numTerms : Int;

    /** The well-separatedness parameter ws. */
    public global val ws : Int;

    // TODO enum - XTENLANG-1118
    public const TIMER_INDEX_TOTAL : Int = 0;
    public const TIMER_INDEX_MULTIPOLE : Int = 1;
    public const TIMER_INDEX_DIRECT : Int = 2;
    public const TIMER_INDEX_COMBINE : Int = 3;
    public const TIMER_INDEX_TRANSFORM : Int = 4;
    public const TIMER_INDEX_FARFIELD : Int = 5;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(6);

    /** All boxes in the octree division of space. */
    private global val boxes : Array[FmmBox](4);

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : Array[ValRail[MMAtom]](1);
    //val atoms : ValRail[MMAtom!];

    /** 
     * A cache of transformations from multipole to local at the same level.
     * Dimensions are:
     * 0: a "virtual" dimension used to replicate the data across all places
     *   // TODO should be done as a global array XTENLANG-787
     * 1: box level in the tree
     * 2: x translation (difference between x coordinate of boxes)
     * 3: y translation
     * 4: z translation
     */
    private global val multipoleTransforms : Array[LocalExpansion](5);

    /** 
     * A cache of multipole translations between parent box centres and child box centres. 
     * Dimensions as per multipoleTransforms
     */
    private global val multipoleTranslations : Array[MultipoleExpansion](5);

    // TODO use shared local variable within getEnergy() - XTENLANG-404
    private var directEnergy : Double = 0.0;
    private var farFieldEnergy : Double = 0.0;

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
                    numAtoms : Int,
                    atoms: Array[ValRail[MMAtom]](1)) {
        val numLevels = Math.max(2, (Math.log(numAtoms / density) / Math.log(8.0) + 1.0 as Int));
        this.numLevels = numLevels;

        var nBox : Int = 0;
        for ((i) in 2..numLevels) {
            nBox += Math.pow2(3*i) as Int;
        }
        val lowestLevelDim = Math.pow2(numLevels);
        this.lowestLevelDim = lowestLevelDim;
        Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + nBox);

        this.numTerms = numTerms;
        this.ws = ws;
        
        this.size = size;

        this.atoms = atoms;

        var boxRegion : Region(4) = [2..2, 0..3, 0..3, 0..3];
        var boxDistribution : Dist(4) = Dist.makeBlock(boxRegion, 0);
        for ((i) in 3..numLevels) {
            val levelDim = Math.pow2(i) as Int;
            val nextLevelRegion : Region(4) = [i..i, 0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val nextLevelDist = Dist.makeBlock(nextLevelRegion, 0);
            boxDistribution = boxDistribution || nextLevelDist;
        }
        Console.OUT.println("boxDistribution: " + boxDistribution);
        
        this.boxes = Array.make[FmmBox](boxDistribution);

        for ((thisLevel) in 2..numLevels-1) {
            //Console.OUT.println("transform level " + thisLevel);
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(4) = [thisLevel..thisLevel, 0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach (boxIndex1 in thisLevelDist) {
                boxes(boxIndex1) = new FmmBox(boxIndex1, numTerms, getParentForChild(boxIndex1));
            }
        }

        
        val lowestLevelRegion : Region(4) = [numLevels..numLevels, 0..lowestLevelDim-1, 0..lowestLevelDim-1, 0..lowestLevelDim-1];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        finish ateach (boxIndex1 in lowestLevelDist) {
            boxes(boxIndex1) = new FmmLeafBox(boxIndex1, numTerms, getParentForChild(boxIndex1));
        }

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
        timer.start(TIMER_INDEX_TOTAL);
        multipoleLowestLevel();
        // direct energy is independent of all subsequent steps of FMM
        val directEnergy = future {getDirectEnergy()};
        combineMultipoles();
        transformToLocal();
        farFieldEnergy = getFarFieldEnergy();
        val totalEnergy = directEnergy() + farFieldEnergy;
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        //Console.OUT.println("multipole lowest level");
        timer.start(TIMER_INDEX_MULTIPOLE);
        //finish ateach (p1 in atoms) {
        finish ateach (p1 in atoms) {
            val localAtoms = atoms(p1);
            foreach ((i) in 0..localAtoms.length-1) {
                val atom = localAtoms(i) as MMAtom!;
                val boxIndex = getLowestLevelBoxIndex(atom);
                async(boxes.dist(boxIndex)) {
                    val remoteAtom = new MMAtom(atom);
                    val leafBox = boxes(boxIndex) as FmmLeafBox!;
                    val boxCentre = leafBox.getCentre(size).sub(remoteAtom.centre);
                    val atomExpansion = MultipoleExpansion.getOlm(remoteAtom.charge, boxCentre, numTerms);
                    atomic {
                        leafBox.atoms.add(remoteAtom);
                        leafBox.multipoleExp.add(atomExpansion);
                    }
                }
            }
        }
        // post-prune leaf boxes
        // TODO prune intermediate empty boxes as well
        val lowestLevelRegion : Region(4) = [numLevels..numLevels, 0..lowestLevelDim-1, 0..lowestLevelDim-1, 0..lowestLevelDim-1];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        finish ateach (boxIndex in lowestLevelDist) {
            val box = boxes(boxIndex) as FmmLeafBox!;
            if (box.atoms.length() == 0) {
                boxes(boxIndex) = null;
            }
        }

        timer.stop(TIMER_INDEX_MULTIPOLE);
    }

    /** 
     * Starting at the bottom level, combines multipole expansions for <= 8 child
     * boxes into a single multipole expansion for the parent box.
     */
    def combineMultipoles() {
        timer.start(TIMER_INDEX_COMBINE);
        for (var level: Int = numLevels; level > 2; level--) {
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val levelDim = Math.pow2(level) as Int;
            val thisLevelRegion : Region(4) = [level..level, 0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach (boxIndex in thisLevelDist) {
                if (boxes(boxIndex) != null) {
                    val child = boxes(boxIndex) as FmmBox!;
                    val childExp = child.multipoleExp;
                    val parent = child.parent;
                    at (parent) {
                        val shift = multipoleTranslations(Point.make([here.id, boxIndex(0), (child.x+1)%2, (child.y+1)%2, (child.z+1)%2])) as MultipoleExpansion!;
                        parent.multipoleExp.translateAndAddMultipole(shift, childExp);
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_COMBINE);
    }

    /**
     * Starting at the top level, for each box, transform multipole
     * expansions for all well-separated boxes (for which the parent
     * box is not also well-separated) to local expansions for this 
     * box.  Translate the local expansion for each parent down to all
     * non-empty child boxes.
     */
    def transformToLocal() {
        //Console.OUT.println("transform level 2");
        timer.start(TIMER_INDEX_TRANSFORM);
        val level2Region : Region(4) = [2..2, 0..3, 0..3, 0..3];
        val level2Dist = boxes.dist | level2Region;
        finish ateach ((level1,x1,y1,z1) in level2Dist) {
            val box1 = boxes(level1,x1,y1,z1) as FmmBox!;
            if (box1 != null) {
                for ((level2,x2,y2,z2) in level2Dist) {
                    if (x1 != x2 || y1 != y2 || z1 != z2) {
                        if (box1.wellSeparated(ws,x2,y2,z2)) {
                            val box2MultipoleExp = getMultipoleExpansionLocalCopy(level2,x2,y2,z2);
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(level2,x2,y2,z2);
                                val transform21 = multipoleTransforms(Point.make([here.id, translation(0), -translation(1), -translation(2), -translation(3)])) as LocalExpansion!;
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                    }
                }
            }
        }
        
        for ((thisLevel) in 3..numLevels) {
            //Console.OUT.println("transform level " + thisLevel);
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(4) = [thisLevel..thisLevel, 0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = boxes.dist | thisLevelRegion;
            finish ateach ((level1,x1,y1,z1) in thisLevelDist) {
                val box1 = boxes(level1,x1,y1,z1) as FmmBox!;
                if (box1 != null) {
                    for ((level2,x2,y2,z2) in thisLevelDist) {
                    if (x1 != x2 || y1 != y2 || z1 != z2) {
                            if (box1.wellSeparated(ws,x2,y2,z2)) {
                                val parentIndex = getParentIndex(level2,x2,y2,z2);
                                if (!box1.parent.wellSeparated(ws, parentIndex)) {
                                    val box2MultipoleExp = getMultipoleExpansionLocalCopy(level2,x2,y2,z2);
                                    if (box2MultipoleExp != null) {
                                        val translation = box1.getTranslationIndex(level2,x2,y2,z2);
                                        val translateP = Point.make([here.id, translation(0), -translation(1), -translation(2), -translation(3)]);
                                        val transform21 = multipoleTransforms(translateP) as LocalExpansion!;
                                        box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                                    }
                                }
                            }
                        }
                    }
                    val box1ParentExp = box1.parent.getLocalExpansionLocalCopy(numTerms);
                    val shift = multipoleTranslations(Point.make([here.id, thisLevel, box1.x%2, box1.y%2, box1.z%2])) as MultipoleExpansion!;
                    box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                }
            }
        }
        timer.stop(TIMER_INDEX_TRANSFORM);
    }

    /**
     * Gets sum of direct (pairwise) energy for all pairs of atoms
     * in non-well-separated boxes. This operations requires only
     * that atoms have already been assigned to boxes, and so can 
     * be done in parallel with other steps of the algorithm.
     */
    def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        // TODO n^2 calculation - to check - remove this
        /*
        var pairwiseEnergy : Double = 0.0;
        for ((i) in 0..(atoms.length - 1)) {
            for ((j) in 0..(i - 1)) {
                val pairEnergy : Double = atoms(j).pairEnergy(atoms(i));
                pairwiseEnergy += 2 * pairEnergy;
            }
        }
        Console.OUT.println("pairwiseEnergy = " + pairwiseEnergy);
        */
        val lowestLevelRegion : Region(4) = [numLevels..numLevels, 0..lowestLevelDim-1, 0..lowestLevelDim-1, 0..lowestLevelDim-1];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        finish ateach (boxIndex1 in lowestLevelDist) {
            val box1 = boxes(boxIndex1) as FmmLeafBox!;
            if (box1 != null) {
                // TODO use shared var - XTENLANG-404
                var thisBoxEnergy : Double = 0.0;
                val length = box1.atoms.length();
                for ((atomIndex1) in 0..length-1) {
                    val atom1 = box1.atoms(atomIndex1);

                    // direct calculation with all atoms in same box
                    for ((sameBoxAtomIndex) in 0..atomIndex1-1) {
                        val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                        val pairEnergy : Double = atom1.pairEnergy(sameBoxAtom);
                        thisBoxEnergy += pairEnergy;
                    }
                }

                // direct in same are counted twice
                thisBoxEnergy *= 2;

                // direct calculation with all atoms in non-well-separated boxes
                for (boxIndex2 in lowestLevelDist) {
                    // TODO halve direct energy! use Morton index
                    if (boxIndex1 != boxIndex2) {
                        if (!box1.wellSeparated(ws, boxIndex2)) {
                            val packedAtoms = at(boxes.dist(boxIndex2)) {getPackedAtomsForBox(boxIndex2)};
                            if (packedAtoms != null) {
                                for (var i : Int = 0; i < packedAtoms.length(); i+=4) {
                                    val atom2Centre = new Point3d(packedAtoms(i+1), packedAtoms(i+2), packedAtoms(i+3));
                                    val atom2Charge = packedAtoms(i);
                                    for ((atomIndex1) in 0..length-1) {
                                        val atom1 = box1.atoms(atomIndex1);
                                        thisBoxEnergy += atom1.charge * atom2Charge / atom1.centre.distance(atom2Centre);
                                    }
                                }
                            }
                        }
                    }
                }
                val thisBoxEnergyFinal = thisBoxEnergy;
                async (this) {atomic {directEnergy += thisBoxEnergyFinal;}}
            }
        }
        timer.stop(TIMER_INDEX_DIRECT);

        return directEnergy;
    }

    /**
     * Calculates the sum of far-field interactions for all atoms:
     * for each atom, calculates the interaction using the local expansion
     * for the enclosing box, which is the combined contribution of
     * all atoms in all well-separated boxes.
     */ 
    def getFarFieldEnergy() : Double {
        //Console.OUT.println("getFarFieldEnergy");
        timer.start(TIMER_INDEX_FARFIELD);
        val lowestLevelRegion : Region(4) = [numLevels..numLevels, 0..lowestLevelDim-1, 0..lowestLevelDim-1, 0..lowestLevelDim-1];
        val lowestLevelDist = boxes.dist | lowestLevelRegion;
        finish ateach (boxIndex1 in lowestLevelDist) {
            val box1 = boxes(boxIndex1) as FmmLeafBox!;
            if (box1 != null) {
                // TODO use shared var - XTENLANG-404
                var thisBoxEnergy : Double = 0.0;
                val length = box1.atoms.length();
                for ((atomIndex1) in 0..length-1) {
                    val atom1 = box1.atoms(atomIndex1);
                    val box1Centre = atom1.centre.sub(box1.getCentre(size));
                    val farFieldEnergy = box1.localExp.getPotential(atom1.charge, box1Centre);
                    thisBoxEnergy += farFieldEnergy;
                }
                val thisBoxEnergyFinal = thisBoxEnergy;
                async (this) {atomic {farFieldEnergy += thisBoxEnergyFinal;}}
            }
        }
        timer.stop(TIMER_INDEX_FARFIELD);

        return farFieldEnergy;
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

    private global def getLowestLevelBoxIndex(atom : MMAtom!) : Point(4) {
        return  Point.make(numLevels, atom.centre.i / size * lowestLevelDim + lowestLevelDim / 2 as Int, atom.centre.j / size * lowestLevelDim + lowestLevelDim / 2 as Int, atom.centre.k / size * lowestLevelDim + lowestLevelDim / 2 as Int);
    }

    private global def getParentForChild(childIndex : Point(4)) : FmmBox {
        if (childIndex(0) == 2)
            // level 2 is highest level
            return null;
        val parentIndex = getParentIndex(childIndex);
        var parent : FmmBox = at (boxes.dist(parentIndex)) {boxes(parentIndex)};
        return parent;
    }

    private global def getParentIndex(boxIndex : Point(4)) : Point(4) {
        return Point.make(boxIndex(0) - 1, boxIndex(1) / 2 , boxIndex(2) / 2, boxIndex(3) / 2);
    }

    private global def getParentIndex(level : Int, x : Int, y : Int, z : Int) : Point(4) {
        return Point.make(level - 1, x / 2 , y / 2, z / 2);
    }

    private global def getPackedAtomsForBox(boxIndex : Point(4)) {
        val box = boxes(boxIndex) as FmmLeafBox!;
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
    private global def getMultipoleExpansionLocalCopy(level : Int, x : Int, y : Int, z : Int) : MultipoleExpansion! {
        val data = at (boxes.dist(level,x,y,z)) {boxes(level,x,y,z) != null? Expansion.getData(numTerms, (boxes(level,x,y,z) as FmmBox!).multipoleExp) : null};
        if (data != null) {
            return new MultipoleExpansion(numTerms, data);
        } else {
            return null;
        }
    }
}

