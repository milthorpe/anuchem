package au.edu.anu.mm;

import x10.util.*;
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
 * For more info see
 * White & Head-Gordon "Derivation and efficient implementation of the fast multipole method", J Chem Phys 101 (8), 1994 
 * and
 * Lashuk et al. "A massively parallel adaptive fast-multipole method on heterogeneous architectures", SC 2009
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
    public const TIMER_INDEX_TREE : Int = 6;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(7);

    /** All boxes in the octree division of space. 
     * ValRail has numLevels-1 elements, for levels [2..numLevels]
     * where boxes(n) = boxes at level n+2
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    private global val boxes : ValRail[Array[FmmBox](3)];

    private global val lowestLevelBoxes : Array[FmmBox](3);

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    private global val atoms : Array[ValRail[MMAtom]](1);

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

    /**
     * An array of locally essential trees (LETs), one for each place.
     */
    private global val locallyEssentialTrees : Array[LocallyEssentialTree](1);

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

        timer.start(TIMER_INDEX_TREE);
        val boxes = constructTree();
        this.boxes = boxes;
        this.lowestLevelBoxes = boxes(numLevels-2);
        this.multipoleTranslations = precomputeTranslations();
        this.multipoleTransforms = precomputeTransforms();
        this.locallyEssentialTrees = createLocallyEssentialTrees();
        timer.stop(TIMER_INDEX_TREE);
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
                async(lowestLevelBoxes.dist(boxIndex)) {
                    val remoteAtom = new MMAtom(atom);
                    val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
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
        finish ateach (boxIndex in lowestLevelBoxes) {
            val box = lowestLevelBoxes(boxIndex) as FmmLeafBox!;
            if (box.atoms.length() == 0) {
                lowestLevelBoxes(boxIndex) = null;
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
            val thisLevel = level;
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelBoxes = boxes(thisLevel-2);
            finish ateach (boxIndex in thisLevelBoxes) {
                if (thisLevelBoxes(boxIndex) != null) {
                    val child = thisLevelBoxes(boxIndex) as FmmBox!;
                    val childExp = child.multipoleExp;
                    val parent = child.parent;
                    at (parent) {
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, (child.x+1)%2, (child.y+1)%2, (child.z+1)%2])) as MultipoleExpansion!;
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

        // start the prefetch of all multipoles required at each place
        prefetchMultipoles();

        val highestLevelBoxes = boxes(0);
        finish ateach (p1 in Dist.makeUnique(highestLevelBoxes.dist.places())) {
            val highestLevelMultipoleCopies = locallyEssentialTrees(here.id).multipoleCopies(0);
            foreach ((x1,y1,z1) in highestLevelBoxes | here) {
                val box1 = highestLevelBoxes(x1,y1,z1) as FmmBox!;
                if (box1 != null) {
                    val vList = box1.getVList();
                    for ((x2,y2,z2) in vList) {
                        // here we force on the multipole value for which a future was previously issued
                        val box2MultipoleExp = highestLevelMultipoleCopies(x2,y2,z2)() as MultipoleExpansion!;
                        if (box2MultipoleExp != null) {
                            val translation = box1.getTranslationIndex(2,x2,y2,z2);
                            val transform21 = multipoleTransforms(Point.make([here.id, translation(0), -translation(1), -translation(2), -translation(3)])) as LocalExpansion!;
                            box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                        }
                    }
                }
            }
        }
        
        for ((thisLevel) in 3..numLevels) {
            //Console.OUT.println("transform level " + thisLevel);
            val thisLevelBoxes = boxes(thisLevel-2);
            finish ateach (p1 in Dist.makeUnique(thisLevelBoxes.dist.places())) {
                val thisLevelMultipoleCopies = locallyEssentialTrees(here.id).multipoleCopies(thisLevel-2);
                foreach ((x1,y1,z1) in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(x1,y1,z1) as FmmBox!;
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for ((x2,y2,z2) in vList) {
                            // force on the multipole value
                            val box2MultipoleExp = thisLevelMultipoleCopies(x2,y2,z2)() as MultipoleExpansion!;
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(thisLevel,x2,y2,z2);
                                val translateP = Point.make([here.id, translation(0), -translation(1), -translation(2), -translation(3)]);
                                val transform21 = multipoleTransforms(translateP) as LocalExpansion!;
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                        val box1ParentExp = box1.parent.getLocalExpansionLocalCopy(numTerms);
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, box1.x%2, box1.y%2, box1.z%2])) as MultipoleExpansion!;
                        box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                    }
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

        // start the prefetch of all atoms required for direct calculations at each place
        prefetchPackedAtoms();
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
        finish ateach (p1 in Dist.makeUnique(lowestLevelBoxes.dist.places())) {
            val packedAtoms = locallyEssentialTrees(here.id).packedAtoms;
            foreach ((x1,y1,z1) in lowestLevelBoxes | here) {
                val box1 = lowestLevelBoxes(x1,y1,z1) as FmmLeafBox!;
                if (box1 != null) {
                    val uList = box1.getUList();

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

                    // direct calculation with all atoms in non-well-separated boxes
                    // interact with "left half" of uList i.e. only boxes with x1<=x1
                    for (box2Index in uList) {
                        if (box2Index(0) < x1 || (box2Index(0) == x1 && box2Index(1) < y1) || (box2Index(0) == x1 && box2Index(1) == y1 && box2Index(2) < z1)) {
                            if (!box1.wellSeparated(ws, box2Index)) {
                                // here we force on the packed atoms for which a future was previously issued
                                val boxAtoms = packedAtoms(box2Index)();
                                if (boxAtoms != null) {
                                    for (var i : Int = 0; i < boxAtoms.length(); i+=4) {
                                        val atom2Centre = new Point3d(boxAtoms(i+1), boxAtoms(i+2), boxAtoms(i+3));
                                        val atom2Charge = boxAtoms(i);
                                        for ((atomIndex1) in 0..length-1) {
                                            val atom1 = box1.atoms(atomIndex1);
                                            thisBoxEnergy += atom1.charge * atom2Charge / atom1.centre.distance(atom2Centre);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    val thisBoxEnergyFinal = 2 * thisBoxEnergy;
                    async (this) {atomic {directEnergy += thisBoxEnergyFinal;}}
                }
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
        finish ateach (boxIndex1 in lowestLevelBoxes) {
            val box1 = lowestLevelBoxes(boxIndex1) as FmmLeafBox!;
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
     * for use in translating multipole expansions from child to parent boxes.
     * This is distributed to all places by replicating the first index in a 
     * cyclic dist across Place.PLACES)
     * TODO workaround due to lack of global immutable arrays - XTENLANG-787
     */
    private def precomputeTranslations() : Array[MultipoleExpansion](5) {
        if (numLevels < 3) {
            return null;
        } else {
            val multipoleTranslations = Array.make[MultipoleExpansion](Dist.makeCyclic([0..Place.MAX_PLACES-1,3..numLevels, 0..1, 0..1, 0..1],0));
            finish ateach ((p) : Point in Dist.makeUnique(Place.places)) {
                for (val(placeId,level,i,j,k) in multipoleTranslations.dist | here) {
                    dim : Int = Math.pow2(level);
                    sideLength : Double = size / dim;
                    val translationVector = new Vector3d((i*2-1) * 0.5 * sideLength,
                                                         (j*2-1) * 0.5 * sideLength,
                                                         (k*2-1) * 0.5 * sideLength);
                    multipoleTranslations(Point.make([placeId, level, i, j, k])) = MultipoleExpansion.getOlm(translationVector as Tuple3d, numTerms);
                }
            }
            return multipoleTranslations;
        }
    }

    /**
     * Precomputes a multipole transform array for use in transforming
     * multipole expansions of well-separated boxes to local expansions
     * at the current box.  This is distributed to all places by replicating 
     * the first index in a cyclic dist across Place.PLACES)
     * TODO workaround due to lack of global immutable arrays - XTENLANG-787
     */
    private def precomputeTransforms() : Array[LocalExpansion](5) {
        var wellSpacedLimit : Region(5) = [0..Place.MAX_PLACES-1,2..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransformRegion : Region(5) = wellSpacedLimit - ([0..Place.MAX_PLACES-1,2..numLevels,-ws..ws,-ws..ws,-ws..ws] as Region);
        //Console.OUT.println("multipoleTransformRegion = " + multipoleTransformRegion);
        val multipoleTransforms = Array.make[LocalExpansion](Dist.makeCyclic(multipoleTransformRegion,0));
        finish ateach ((p) : Point in Dist.makeUnique(Place.places)) {
            for (val(placeId,level,i,j,k) in multipoleTransforms.dist | here) {
                dim : Int = Math.pow2(level);
                sideLength : Double = size / dim;
                val translationVector = new Vector3d(i * sideLength,
                                                 j * sideLength,
                                                 k * sideLength);
                multipoleTransforms(Point.make([placeId, level, i, j, k])) = LocalExpansion.getMlm(translationVector, numTerms);
            }
        }
        return multipoleTransforms;
    }

    private def constructTree() : ValRail[Array[FmmBox](3)] {
        val boxesTemp = Rail.make[Array[FmmBox](3)](numLevels-1); // there's no level 1
        for ((thisLevel) in 2..numLevels) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(3) = [0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = Dist.makeBlock(thisLevelRegion, 0);
            boxesTemp(thisLevel-2) = Array.make[FmmBox](thisLevelDist);
            Console.OUT.println("level " + thisLevel + " dist: " + thisLevelDist);
        }
        val boxesValRail = boxesTemp as ValRail[Array[FmmBox](3)];

        for ((thisLevel) in 2..numLevels-1) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelBoxes = boxesValRail(thisLevel-2);
            finish ateach ((x,y,z) in thisLevelBoxes) {
                val box = new FmmBox(thisLevel, x, y, z, numTerms, getParentForChild(boxesValRail, thisLevel, x,y,z));
                createVList(box);
                thisLevelBoxes(x,y,z) = box;
            }
        }

        val lowestLevelBoxes = boxesValRail(numLevels-2);
        finish ateach ((x,y,z) in lowestLevelBoxes) {
            val box = new FmmLeafBox(numLevels, x, y, z, numTerms, getParentForChild(boxesValRail, numLevels, x,y,z));
            createUList(box);
            createVList(box);
            lowestLevelBoxes(x,y,z) = box;
        }

        return boxesValRail;
    }

    /**
     * Creates the locally essential tree at each place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private def createLocallyEssentialTrees() : Array[LocallyEssentialTree](1) {
        val locallyEssentialTrees = Array.make[LocallyEssentialTree](Dist.makeUnique(), (Point)=> null);
        finish ateach ((p1) in locallyEssentialTrees) {
            val uMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
            val uMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
            val combinedUSet = new HashSet[Point(3)]();
            for (boxIndex1 in lowestLevelBoxes | here) {
                val box1 = lowestLevelBoxes(boxIndex1) as FmmLeafBox!;
                if (box1 != null) {
                    val uList = box1.getUList();
                    for (boxIndex2 in uList) {
                        if (combinedUSet.add(boxIndex2)) {
                            for ((i) in 0..2) {
                                uMin(i) = Math.min(uMin(i), boxIndex2(i));
                                uMax(i) = Math.max(uMax(i), boxIndex2(i));        
                            }
                        }
                    }
                }
            }
            val uListMin = uMin as ValRail[Int](3);
            val uListMax = uMax as ValRail[Int](3);
            val combinedUList = combinedUSet.toValRail();

            val combinedVList = Rail.make[ValRail[Point(3)]](numLevels-1);
            val vListMin = Rail.make[ValRail[Int](3)](numLevels-1);
            val vListMax = Rail.make[ValRail[Int](3)](numLevels-1);
            for ((thisLevel) in 2..numLevels) {
                val vMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
                val vMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
                //Console.OUT.println("create combined V-list for level " + thisLevel + " at " + here);
                val combinedVSet = new HashSet[Point(3)]();
                val thisLevelBoxes = boxes(thisLevel-2);
                for (boxIndex1 in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(boxIndex1) as FmmBox!;
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for (boxIndex2 in vList) {
                            if (combinedVSet.add(boxIndex2)) {
                                for ((i) in 0..2) {
                                    vMin(i) = Math.min(vMin(i), boxIndex2(i));
                                    vMax(i) = Math.max(vMax(i), boxIndex2(i));        
                                }
                            }
                        }
                    }
                }
                //Console.OUT.println("done " + combinedVSet.size());
                vListMin(thisLevel-2) = vMin as ValRail[Int](3);
                vListMax(thisLevel-2) = vMax as ValRail[Int](3);
                combinedVList(thisLevel-2) = combinedVSet.toValRail();
            }
            locallyEssentialTrees(p1) = new LocallyEssentialTree(combinedUList,
                                                                 combinedVList as ValRail[ValRail[Point(3)]],
                                                                 uListMin,
                                                                 uListMax,
                                                                 vListMin as ValRail[ValRail[Int]],
                                                                 vListMax as ValRail[ValRail[Int]]);
        }
        return locallyEssentialTrees;
    }

    /**
     * "Pre-fetches" the multipole copies for the V-List of the locally
     * essential tree at each place, using futures.
     */
    private def prefetchMultipoles() {
        finish ateach (p1 in locallyEssentialTrees) {
            val myLET = locallyEssentialTrees(p1) as LocallyEssentialTree!;
            val combinedVList = myLET.combinedVList;
            for ((level) in 2..numLevels) {
                val thisLevelBoxes = boxes(level-2);
                for ((x,y,z) in combinedVList(level-2)) {
                    myLET.multipoleCopies(level-2)(x,y,z) = future {getMultipoleExpansionLocalCopy(thisLevelBoxes,x,y,z)};
                }
            }
        }
    }

    /**
     * "Pre-fetches" the packed atoms for the U-List of the locally
     * essential tree at each place, using futures.
     */
    private def prefetchPackedAtoms() {
        finish ateach (p1 in locallyEssentialTrees) {
            val myLET = locallyEssentialTrees(p1) as LocallyEssentialTree!;
            val myCombinedUList = myLET.combinedUList;
            for ((x,y,z) in myCombinedUList) {
                myLET.packedAtoms(x,y,z) = future {at (lowestLevelBoxes.dist(x,y,z)) {getPackedAtomsForBox(x, y, z)}};
            }
        }
    }

    private global def getLowestLevelBoxIndex(atom : MMAtom!) : Point(3) {
        return  Point.make(atom.centre.i / size * lowestLevelDim + lowestLevelDim / 2 as Int, atom.centre.j / size * lowestLevelDim + lowestLevelDim / 2 as Int, atom.centre.k / size * lowestLevelDim + lowestLevelDim / 2 as Int);
    }

    private global def getParentForChild(boxes : ValRail[Array[FmmBox](3)], level : Int, x : Int, y : Int, z : Int) : FmmBox {
        if (level == 2)
            // level 2 is highest level
            return null;
        val parentIndex = getParentIndex(x, y, z);
        var parent : FmmBox = at (boxes(level-3).dist(parentIndex)) {boxes(level-3)(parentIndex)};
        return parent;
    }

    private global def getParentIndex(x : Int, y : Int, z : Int) : Point(3) {
        return Point.make(x / 2 , y / 2, z / 2);
    }

    private global def getPackedAtomsForBox(x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox!;
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
    private global def getMultipoleExpansionLocalCopy(thisLevelBoxes : Array[FmmBox](3), x : Int, y : Int, z : Int) : MultipoleExpansion! {
        val data = at (thisLevelBoxes.dist(x,y,z)) {thisLevelBoxes(x,y,z) != null? Expansion.getData(numTerms, (thisLevelBoxes(x,y,z) as FmmBox!).multipoleExp) : null};
        if (data != null) {
            return new MultipoleExpansion(numTerms, data);
        } else {
            return null;
        }
    }

    /**
     * Creates the U-list of <code>box</code>.
     * The U-list consists of all leaf boxes not well-separated from <code>box</code>.
     */
    private global def createUList(box : FmmLeafBox!) {
        val uList = new GrowableRail[Point(3)]();
        for ((x) in Math.max(0,box.x-ws)..Math.min(lowestLevelDim-1,box.x+ws)) {
            for ((y) in Math.max(0,box.y-ws)..Math.min(lowestLevelDim-1,box.y+ws)) {
                for ((z) in Math.max(0,box.z-ws)..Math.min(lowestLevelDim-1,box.z+ws)) {
                    uList.add(Point.make(x,y,z));
                }
            }
        }
        box.setUList(uList.toValRail());
    }

    /**
     * Creates the V-list of <code>box</code>.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent of <code>box</code>.
     */
    private global def createVList(box : FmmBox!) {
        val levelDim = Math.pow2(box.level);
        val xOffset = box.x%2 == 1 ? -1 : 0;
        val yOffset = box.y%2 == 1 ? -1 : 0;
        val zOffset = box.z%2 == 1 ? -1 : 0;
        val vList = new GrowableRail[Point(3)]();
        for ((x) in Math.max(0,box.x-2*ws+xOffset)..Math.min(levelDim-1,box.x+2*ws+1+xOffset)) {
            for ((y) in Math.max(0,box.y-2*ws+yOffset)..Math.min(levelDim-1,box.y+2*ws+1+yOffset)) {
                for ((z) in Math.max(0,box.z-2*ws+zOffset)..Math.min(levelDim-1,box.z+2*ws+1+zOffset)) {
                    if (box.wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setVList(vList.toValRail());
    }
}

