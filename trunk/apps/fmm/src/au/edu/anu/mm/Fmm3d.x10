/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.util.*;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;

/**
 * This class implements the Fast Multipole Method for electrostatic
 * calculations in a cubic simulation space centred at the origin.
 * <p>
 * The 3D simulation space is divided in an octree of <code>numLevels</code> levels.
 * </p> 
 * @see White & Head-Gordon (1994). "Derivation and efficient implementation 
 *  of the fast multipole method". J Chem Phys 101 (8)
 * @see Lashuk et al. (2009). "A massively parallel adaptive fast-multipole 
 *  method on heterogeneous architectures". Proceedings of SC 2009
 * @author milthorpe
 */
public class Fmm3d {
    /** The number of levels in the octree. */
    public val numLevels : Int;

    /** The number of lowest level boxes along one side of the cube. */
    public val lowestLevelDim : Int;

    /** 
     * To ensure balanced rounding errors within the multipole and local 
     * calculations, all force/energy calculations are performed within 
     * an offset cube with top-left-front corner (-size/2, -size/2, -size/2).
     * This is the offset vector from the "real" (input) cube top-left-front
     * corner to the FMM top-left-front corner. 
     */
    public val offset : Vector3d;

    /** The side length of the cube. */
    public val size : Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public val numTerms : Int;

    /** The well-separatedness parameter ws. */
    public val ws : Int;

    /** 
     * Return the top level of boxes actually used in the method.
     * This is 0 for the periodic FMM and 2 for the non-periodic FMM.
     */
    protected val topLevel : Int;

    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int = 0;
    public static val TIMER_INDEX_PREFETCH : Int = 1;
    public static val TIMER_INDEX_DIRECT : Int = 2;
    public static val TIMER_INDEX_MULTIPOLE : Int = 3;
    public static val TIMER_INDEX_COMBINE : Int = 4;
    public static val TIMER_INDEX_TRANSFORM : Int = 5;
    public static val TIMER_INDEX_FARFIELD : Int = 6;
    public static val TIMER_INDEX_TREE : Int = 7;
    public static val TIMER_INDEX_SETUP : Int = 8;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(9);

    /** All boxes in the octree division of space. 
     * ValRail has numLevels elements, for levels [1..numLevels]
     * where boxes(n) = boxes at level n+1
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    val boxes : ValRail[DistArray[FmmBox](3){rect}];

    val lowestLevelBoxes : DistArray[FmmBox](3){rect};

    /** The atoms in the simulation, divided up into an array of ValRails, one for each place. */
    protected val atoms : DistArray[ValRail[MMAtom]](1);

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
    val multipoleTransforms : DistArray[LocalExpansion](5){rect};

    /** 
     * A cache of multipole translations between parent box centres and child box centres. 
     * Dimensions as per multipoleTransforms
     */
    val multipoleTranslations : DistArray[MultipoleExpansion](5){rect};

    /**
     * An array of locally essential trees (LETs), one for each place.
     */
    protected val locallyEssentialTrees : DistArray[LocallyEssentialTree](1);

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
                    topLeftFront : Point3d,
                    size : Double,  
                    numAtoms : Int,
                    atoms: DistArray[ValRail[MMAtom]](1)) {
        // topLevel in regular FMM is 2 (boxes higher than this cannot be well-spaced)
        this(density, numTerms, ws, topLeftFront, size, numAtoms, atoms, 2);
    }

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     * @param topLevel the topmost level for which boxes in the octree are used
     */
    protected def this(density : Double, 
                    numTerms : Int,
                    ws : Int,
                    topLeftFront : Point3d,
                    size : Double,  
                    numAtoms : Int,
                    atoms: DistArray[ValRail[MMAtom]](1),
                    topLevel : Int) {
        this.topLevel = topLevel;
        val numLevels = Math.max(2, (Math.log(numAtoms / density) / Math.log(8.0) + 1.0) as Int);
        this.numLevels = numLevels;

        var nBox : Int = 0;
        for ([i] in topLevel..numLevels) {
            nBox += Math.pow2(3*i) as Int;
        }
        val lowestLevelDim = Math.pow2(numLevels);
        this.lowestLevelDim = lowestLevelDim;
        Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + nBox);

        this.numTerms = numTerms;
        this.ws = ws;

        val offset = Point3d(-size/2.0, -size/2.0, -size/2.0) - topLeftFront;
        this.offset = offset;
        this.size = size;

        this.atoms = atoms;

        timer.start(TIMER_INDEX_TREE);
        val boxes = constructTree();
        this.boxes = boxes;
        this.lowestLevelBoxes = boxes(numLevels);
        this.multipoleTranslations = precomputeTranslations();
        this.multipoleTransforms = precomputeTransforms();
        this.locallyEssentialTrees = createLocallyEssentialTrees();
        assignAtomsToBoxes(atoms, boxes(numLevels));
        timer.stop(TIMER_INDEX_TREE);
    }
    
    public def calculateEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);
        finish {
            async {
                prefetchPackedAtoms();
            }
            multipoleLowestLevel();
            combineMultipoles();
            transformToLocal();
        }
        val totalEnergy = getDirectEnergy() + getFarFieldEnergy();
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }


    private def assignAtomsToBoxes(atoms: DistArray[ValRail[MMAtom]](1), lowestLevelBoxes : DistArray[FmmBox](3)) {
        //Console.OUT.println("assignAtomsToBoxes");
        finish ateach (p1 in atoms) {
            val localAtoms = atoms(p1);
            finish foreach ([i] in 0..localAtoms.length-1) {
                val atom = localAtoms(i);
                val charge = atom.charge;
                val offsetCentre = atom.centre + offset;
                val boxIndex = getLowestLevelBoxIndex(offsetCentre);
                at(lowestLevelBoxes.dist(boxIndex)) {
                    val remoteAtom = new MMAtom(offsetCentre, charge);
                    val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox;
                    atomic {
                        leafBox.atoms.add(remoteAtom);
                    }
                }
            }
        }
        // post-prune leaf boxes
        // TODO prune intermediate empty boxes as well
        finish ateach (boxIndex in lowestLevelBoxes) {
            val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
            if (box.atoms.length() == 0) {
                lowestLevelBoxes(boxIndex) = null;
            }
        }
    }

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        //Console.OUT.println("multipole lowest level");
        timer.start(TIMER_INDEX_MULTIPOLE);
        finish ateach (boxIndex in lowestLevelBoxes) {
            val leafBox = lowestLevelBoxes(boxIndex) as FmmLeafBox;
            if (leafBox != null) {
                val boxLocation = leafBox.getCentre(size);
                for ([i] in 0..leafBox.atoms.length()-1) {
                    val atom = leafBox.atoms(i);
                    val atomLocation = leafBox.getCentre(size).vector(atom.centre);
                    val atomExpansion = MultipoleExpansion.getOlm(atom.charge, atomLocation, numTerms);
                    leafBox.multipoleExp.add(atomExpansion);
                }
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
        for (var level: Int = numLevels; level > topLevel; level--) {
            val thisLevel = level;
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelBoxes = boxes(thisLevel);
            finish ateach (boxIndex in thisLevelBoxes) {
                val child = thisLevelBoxes(boxIndex);
                if (child != null) {
                    val childExp = child.multipoleExp;
                    val parent = child.parent;
                    at (parent) {
                        val shift = multipoleTranslations(Point.make([here.id, thisLevel, (child.x+1)%2, (child.y+1)%2, (child.z+1)%2]));
                        parent().multipoleExp.translateAndAddMultipole(shift, childExp);
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
        timer.start(TIMER_INDEX_TRANSFORM);
        
        for ([thisLevel] in (topLevel)..numLevels) {
            //Console.OUT.println("transform level " + thisLevel);
            val thisLevelBoxes = boxes(thisLevel);
            finish ateach (p1 in Dist.makeUnique(thisLevelBoxes.dist.places())) {
                val myLET = locallyEssentialTrees(p1);
                val combinedVList = myLET.combinedVList;

                val thisLevelMultipoleCopies = myLET.multipoleCopies(thisLevel);
                if (thisLevel == topLevel) {
                    // must fetch top level multipoles synchronously before starting
                    val thisLevelBoxes = boxes(thisLevel);
                    finish foreach ([x,y,z] in combinedVList(thisLevel)) {
                        thisLevelMultipoleCopies(x,y,z) = getMultipoleExpansionLocalCopy(thisLevelBoxes,x,y,z);
                    }
                }

                // can fetch next level multipoles asynchronously while computing this level
                async {
                    val lowerLevel = thisLevel+1;
                    if (lowerLevel <= numLevels) {
                        val lowerLevelCopies = myLET.multipoleCopies(lowerLevel);
                        val lowerLevelBoxes = boxes(lowerLevel);
                        for ([x,y,z] in combinedVList(lowerLevel)) {
                            lowerLevelCopies(x,y,z) = getMultipoleExpansionLocalCopy(lowerLevelBoxes,x,y,z);
                        }
                    }
                }

                finish foreach ([x1,y1,z1] in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(x1,y1,z1);
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for ([x2,y2,z2] in vList) {
                            // force on the multipole value
                            val box2MultipoleExp = thisLevelMultipoleCopies(x2,y2,z2);
                            if (box2MultipoleExp != null) {
                                val translation = box1.getTranslationIndex(thisLevel,x2,y2,z2);
                                val translateP = Point.make([here.id, translation(0), translation(1), translation(2), translation(3)]);
                                val transform21 = multipoleTransforms(translateP);
                                box1.localExp.transformAndAddToLocal(transform21, box2MultipoleExp);
                            }
                        }
                        if (thisLevel > topLevel) {
                            val box1Parent = box1.parent;
                            val box1ParentExp = at (box1Parent) {box1Parent().localExp};
                            val shift = multipoleTranslations(Point.make([here.id, thisLevel, box1.x%2, box1.y%2, box1.z%2]));
                            box1.localExp.translateAndAddLocal(shift, box1ParentExp);
                        }
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_TRANSFORM);
    }

    /**
     * Fetch all atoms required for direct calculations at each place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchPackedAtoms() : void {
        timer.start(TIMER_INDEX_PREFETCH);
        finish ateach (p1 in locallyEssentialTrees) {
            val myLET = locallyEssentialTrees(p1);
            val myCombinedUList = myLET.combinedUList;
            val packedAtoms = myLET.packedAtoms;

            val uListPlaces = new HashMap[Int,GrowableRail[Point(3)]](26); // a place may have up to 26 immediate neighbours
            
            // separate the uList into partial lists stored at each nearby place
            for (boxIndex in myCombinedUList) {
                val placeId = lowestLevelBoxes.dist(boxIndex).id;
                var uListForPlace : GrowableRail[Point(3)] = uListPlaces.getOrElse(placeId, null);
                if (uListForPlace == null) {
                    uListForPlace = new GrowableRail[Point(3)]();
                    uListPlaces.put(placeId, uListForPlace);
                }
                uListForPlace.add(boxIndex);
            }

            // retrieve the partial list for each place and store into my LET
            finish foreach(placeEntry in uListPlaces.entries()) {
                val placeId = placeEntry.getKey();
                val uListForPlace = placeEntry.getValue();
                val uListValRail = uListForPlace.toValRail();
                val packedForPlace = at (Place.place(placeId)) { getPackedAtomsForBoxList(uListValRail)};
                for ([i] in 0..uListValRail.length()-1) {
                    myLET.packedAtoms(uListValRail(i)) = packedForPlace(i);
                }
            }
        }
        timer.stop(TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of box indices as Point(3) stored at a single
     * place, returns a ValRail, each element of which is in turn
     * a ValRail of MMAtom.PackedRepresentation containing the 
     * packed atoms for each box.
     */
    private def getPackedAtomsForBoxList(boxList : ValRail[Point(3)]) {
        val packedAtomList = ValRail.make[ValRail[MMAtom.PackedRepresentation]](boxList.length(), 
                                                            (i : Int) => 
                                                                getPackedAtomsForBox(boxList(i)(0), 
                                                                                     boxList(i)(1),
                                                                                     boxList(i)(2))
                                                            );
        return packedAtomList;
    }

    /**
     * Gets sum of direct (pairwise) energy for all pairs of atoms
     * in non-well-separated boxes. This operations requires only
     * that atoms have already been assigned to boxes, and so can 
     * be done in parallel with other steps of the algorithm.
     */
    def getDirectEnergy() : Double {
        //Console.OUT.println("direct");
        timer.start(TIMER_INDEX_DIRECT);

        val directEnergy = finish (SumReducer()) {
            ateach (p1 in locallyEssentialTrees) {
                val myLET = locallyEssentialTrees(p1);
                val packedAtoms = myLET.packedAtoms;
                var thisPlaceEnergy : Double = 0.0;
                for ([x1,y1,z1] in lowestLevelBoxes | here) {
                    val box1 = lowestLevelBoxes(x1,y1,z1) as FmmLeafBox;
                    if (box1 != null) {
                        val length = box1.atoms.length();
                        for ([atomIndex1] in 0..length-1) {
                            val atom1 = box1.atoms(atomIndex1);

                            // direct calculation with all atoms in same box
                            for ([sameBoxAtomIndex] in 0..atomIndex1-1) {
                                val sameBoxAtom = box1.atoms(sameBoxAtomIndex);
                                val pairEnergy = atom1.charge * sameBoxAtom.charge / atom1.centre.distance(sameBoxAtom.centre);
                                thisPlaceEnergy += pairEnergy;
                            }

                            // direct calculation with all atoms in non-well-separated boxes
                            val uList = box1.getUList();
                            for ([x2,y2,z2] in uList) {
                                val boxAtoms = packedAtoms(x2,y2,z2);
                                if (boxAtoms != null) {
                                    for ([otherBoxAtomIndex] in 0..(boxAtoms.length()-1)) {
                                        val atom2Packed = boxAtoms(otherBoxAtomIndex);
                                        thisPlaceEnergy += atom1.charge * atom2Packed.charge / atom1.centre.distance(atom2Packed.centre);
                                    }
                                }
                            }
                        }
                    }
                }
                offer thisPlaceEnergy;
            }
        };
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
        val farFieldEnergy = finish(SumReducer()) {
            ateach (boxIndex1 in lowestLevelBoxes) {
                val box1 = lowestLevelBoxes(boxIndex1) as FmmLeafBox;
                if (box1 != null) {
                    var thisBoxEnergy : Double = 0.0;
                    val length = box1.atoms.length();
                    for ([atomIndex1] in 0..length-1) {
                        val atom1 = box1.atoms(atomIndex1);
                        val box1Centre = atom1.centre.vector(box1.getCentre(size));
                        thisBoxEnergy += box1.localExp.getPotential(atom1.charge, box1Centre);
                    }
                    offer thisBoxEnergy;
                }
            }
        };
        timer.stop(TIMER_INDEX_FARFIELD);

        return farFieldEnergy / 2.0;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public def apply(a:Double, b:Double) = (a + b);
    }

    /**
     * Precomputes multipole translations and multipole-to-local transformations,
     * for use in translating multipole expansions from child to parent boxes.
     * This is distributed to all places by replicating the first index in a 
     * block dist across Place.PLACES)
     * TODO workaround due to lack of global immutable arrays - XTENLANG-787
     */
    private def precomputeTranslations() : DistArray[MultipoleExpansion](5){rect} {
        val topChildLevel = topLevel + 1;
        if (numLevels < topChildLevel) {
            return null;
        } else {
            val multipoleTranslations = DistArray.make[MultipoleExpansion](Dist.makeBlock([0..Place.MAX_PLACES-1,topChildLevel..numLevels, 0..1, 0..1, 0..1],0));
            finish ateach ([placeId,level,i,j,k] in multipoleTranslations) {
                val dim = Math.pow2(level);
                val sideLength = size / dim;
                val translationVector = Vector3d((i*2-1) * 0.5 * sideLength,
                                                 (j*2-1) * 0.5 * sideLength,
                                                 (k*2-1) * 0.5 * sideLength);
                multipoleTranslations(Point.make([placeId, level, i, j, k])) = MultipoleExpansion.getOlm(translationVector, numTerms);
            }
            return multipoleTranslations;
        }
    }

    /**
     * Precomputes a multipole transform array for use in transforming
     * multipole expansions of well-separated boxes to local expansions
     * at the current box.  This is distributed to all places by replicating 
     * the first index in a block dist across Place.PLACES)
     * TODO workaround due to lack of global immutable arrays - XTENLANG-787
     */
    private def precomputeTransforms() : DistArray[LocalExpansion](5){rect} {
        var multipoleTransformRegion : Region(5){rect} = [0..Place.MAX_PLACES-1,topLevel..numLevels,-(ws+3)..ws+3,-(ws+3)..ws+3,-(ws+3)..ws+3];
        val multipoleTransforms = DistArray.make[LocalExpansion](Dist.makeBlock(multipoleTransformRegion,0));
        finish ateach ([placeId,level,i,j,k] in multipoleTransforms) {
            val dim = Math.pow2(level);
            val sideLength = size / dim;
            val translationVector = Vector3d(i * sideLength,
                                             j * sideLength,
                                             k * sideLength);
            multipoleTransforms(Point.make([placeId, level, i, j, k])) = LocalExpansion.getMlm(translationVector, numTerms);
        }
        return multipoleTransforms;
    }

    private def constructTree() : ValRail[DistArray[FmmBox](3){rect}] {
        val boxesTemp = Rail.make[DistArray[FmmBox](3){rect}](numLevels+1);
        for ([thisLevel] in topLevel..numLevels) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelRegion : Region(3){rect} = [0..levelDim-1, 0..levelDim-1, 0..levelDim-1];
            val thisLevelDist = MortonDist.make(thisLevelRegion);
            boxesTemp(thisLevel) = DistArray.make[FmmBox](thisLevelDist);
            Console.OUT.println("level " + thisLevel + " dist: " + thisLevelDist);
        }
        val boxesValRail = ValRail.make(boxesTemp);

        for ([thisLevel] in topLevel..numLevels-1) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelBoxes = boxesValRail(thisLevel);
            finish ateach ([x,y,z] in thisLevelBoxes) {
                val box = new FmmBox(thisLevel, x, y, z, numTerms, getParentForChild(boxesValRail, thisLevel, x,y,z));
                createVList(box);
                thisLevelBoxes(x,y,z) = box;
            }
        }

        val lowestLevelBoxes = boxesValRail(numLevels);
        finish ateach ([x,y,z] in lowestLevelBoxes) {
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
    private def createLocallyEssentialTrees() : DistArray[LocallyEssentialTree](1) {
        val locallyEssentialTrees = DistArray.make[LocallyEssentialTree](Dist.makeUnique(), (Point)=> null);
        finish ateach ([p1] in locallyEssentialTrees) {
            val uMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
            val uMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
            val combinedUSet = new HashSet[Point(3)]();
            for ([x,y,z] in lowestLevelBoxes | here) {
                val box1 = lowestLevelBoxes(x,y,z) as FmmLeafBox;
                if (box1 != null) {
                    val uList = box1.getUList();
                    for (boxIndex2 in uList) {
                        if (combinedUSet.add(boxIndex2)) {
                            for ([i] in 0..2) {
                                uMin(i) = Math.min(uMin(i), boxIndex2(i));
                                uMax(i) = Math.max(uMax(i), boxIndex2(i));        
                            }
                        }
                    }
                }
            }
            val uListMin = ValRail.make(uMin);
            val uListMax = ValRail.make(uMax);
            val combinedUList = combinedUSet.toValRail();

            val combinedVList = Rail.make[ValRail[Point(3)]](numLevels+1);
            val vListMin = Rail.make[ValRail[Int](3)](numLevels+1);
            val vListMax = Rail.make[ValRail[Int](3)](numLevels+1);
            for ([thisLevel] in topLevel..numLevels) {
                val vMin = Rail.make[Int](3, (Int) => Int.MAX_VALUE);
                val vMax = Rail.make[Int](3, (Int) => Int.MIN_VALUE);
                //Console.OUT.println("create combined V-list for level " + thisLevel + " at " + here);
                val combinedVSet = new HashSet[Point(3)]();
                val thisLevelBoxes = boxes(thisLevel);
                for ([x,y,z] in thisLevelBoxes | here) {
                    val box1 = thisLevelBoxes(x,y,z);
                    if (box1 != null) {
                        val vList = box1.getVList();
                        for (boxIndex2 in vList) {
                            if (combinedVSet.add(boxIndex2)) {
                                for ([i] in 0..2) {
                                    vMin(i) = Math.min(vMin(i), boxIndex2(i));
                                    vMax(i) = Math.max(vMax(i), boxIndex2(i));        
                                }
                            }
                        }
                    }
                }
                //Console.OUT.println("done " + combinedVSet.size());
                vListMin(thisLevel) = ValRail.make(vMin);
                vListMax(thisLevel) = ValRail.make(vMax);
                combinedVList(thisLevel) = combinedVSet.toValRail();
            }
            locallyEssentialTrees(p1) = new LocallyEssentialTree(combinedUList,
                                                                 ValRail.make(combinedVList),
                                                                 uListMin,
                                                                 uListMax,
                                                                 ValRail.make(vListMin),
                                                                 ValRail.make(vListMax));
        }
        return locallyEssentialTrees;
    }

    private def getLowestLevelBoxIndex(offsetCentre : Point3d) : Point(3) {
        return  Point.make((offsetCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as Int);
    }

    private def getParentForChild(boxes : ValRail[DistArray[FmmBox](3)], level : Int, x : Int, y : Int, z : Int) : GlobalRef[FmmBox] {
        if (level == topLevel)
            return GlobalRef[FmmBox](null);
        return (at (boxes(level-1).dist(x/2, y/2, z/2)) {GlobalRef[FmmBox](boxes(level-1)(x/2, y/2, z/2))});
    }

    private def getPackedAtomsForBox(x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox;
        if (box != null) {
            return box.getPackedAtoms();
        } else {
            return null;
        }
    }

    /**
     * @return a local copy at the current place of a box's multipole expansion
     */
    private def getMultipoleExpansionLocalCopy(thisLevelBoxes : DistArray[FmmBox](3), x : Int, y : Int, z : Int) : MultipoleExpansion {
        return at (thisLevelBoxes.dist(x,y,z)) {thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null};
    }

    /**
     * Creates the U-list of <code>box</code>.
     * The U-list consists of all leaf boxes not well-separated from <code>box</code>.
     */
    private def createUList(box : FmmLeafBox) {
        // interact with "left half" of uList i.e. only boxes with x<=box.x
        val uList = new GrowableRail[Point(3)]();
        for ([x] in Math.max(0,box.x-ws)..box.x) {
            for ([y] in Math.max(0,box.y-ws)..Math.min(lowestLevelDim-1,box.y+ws)) {
                for ([z] in Math.max(0,box.z-ws)..Math.min(lowestLevelDim-1,box.z+ws)) {
                    if (x < box.x || (x == box.x && y < box.y) || (x == box.x && y == box.y && z < box.z)) {
                        uList.add(Point.make(x,y,z));
                    }
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
    private def createVList(box : FmmBox) {
        val levelDim = Math.pow2(box.level);
        val xOffset = box.x%2 == 1 ? -1 : 0;
        val yOffset = box.y%2 == 1 ? -1 : 0;
        val zOffset = box.z%2 == 1 ? -1 : 0;
        val vList = new GrowableRail[Point(3)]();
        for ([x] in Math.max(0,box.x-2*ws+xOffset)..Math.min(levelDim-1,box.x+2*ws+1+xOffset)) {
            for ([y] in Math.max(0,box.y-2*ws+yOffset)..Math.min(levelDim-1,box.y+2*ws+1+yOffset)) {
                for ([z] in Math.max(0,box.z-2*ws+zOffset)..Math.min(levelDim-1,box.z+2*ws+1+zOffset)) {
                    if (box.wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        box.setVList(vList.toValRail());
    }
}

