/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2011.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.HashSet;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.PointCharge;
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
    public static val TIMER_INDEX_DOWNWARD : Int = 5;
    public static val TIMER_INDEX_TREE : Int = 6;
    public static val TIMER_INDEX_PLACEHOLDER : Int = 7;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(8);

    /** All boxes in the octree division of space. 
     * Array has numLevels elements, for levels [1..numLevels]
     * where boxes(n) = boxes at level n+1
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    val boxes : Rail[DistArray[FmmBox](3)];

    /** 
     * The locally essential tree at each place. 
     * @see Lashuk et al. (2009).
     */
    protected val locallyEssentialTree : PlaceLocalHandle[LocallyEssentialTree];

    /**
     * Are boundary conditions periodic?
     */
    val periodic : boolean;

    /** A local cache at every place of the operator arrays for transformations and translations. */
    val fmmOperators : PlaceLocalHandle[FmmOperators];

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
                    numAtoms : Int) {
        // topLevel in regular FMM is 2 (boxes higher than this cannot be well-spaced)
        this(density, numTerms, ws, topLeftFront, size, numAtoms, 2, false);
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
                    topLevel : Int,
                    periodic : boolean) {
        this.topLevel = topLevel;
        this.periodic = periodic;
        val numLevels = Math.max(2, (Math.log(numAtoms / density) / Math.log(8.0) + 1.0) as Int);
        this.numLevels = numLevels;

        var nBox : Int = 0;
        for (i in topLevel..numLevels) {
            nBox += Math.pow2(3*i) as Int;
        }
        val lowestLevelDim = Math.pow2(numLevels);
        this.lowestLevelDim = lowestLevelDim;
        //Console.OUT.println("numLevels = " + numLevels + " maxBoxes = " + nBox);

        this.numTerms = numTerms;
        this.ws = ws;

        val offset = Point3d(-size/2.0, -size/2.0, -size/2.0) - topLeftFront;
        this.offset = offset;
        this.size = size;

        timer.start(TIMER_INDEX_TREE);
        val boxes = constructTree(numLevels, topLevel, numTerms, ws, periodic);
        this.boxes = boxes;
        this.locallyEssentialTree = PlaceLocalHandle.make[LocallyEssentialTree](Dist.makeUnique(), () => Fmm3d.createLocallyEssentialTree(numLevels, topLevel, boxes, periodic));

        this.fmmOperators = PlaceLocalHandle.make[FmmOperators](Dist.makeUnique(), () => new FmmOperators(numTerms, ws));

        timer.stop(TIMER_INDEX_TREE);
    }
    
    public def calculateEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);
        val farFieldEnergy:Double;
        finish {
            async {
                prefetchRemoteAtoms();
            }
            upwardPass();
            farFieldEnergy = downwardPass();
        }
        val totalEnergy = getDirectEnergy() + farFieldEnergy;
        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }


    public def assignAtomsToBoxes(atoms:DistArray[Rail[MMAtom]](1)) {
        timer.start(TIMER_INDEX_TREE);
        //Console.OUT.println("assignAtomsToBoxes");
        val lowestLevelBoxes = boxes(numLevels);
        val offset = this.offset; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelDim = this.lowestLevelDim; // TODO shouldn't be necessary XTENLANG-1913
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val boxAtomsTemp = DistArray.make[ArrayList[PointCharge]](lowestLevelBoxes.dist, (Point) => new ArrayList[PointCharge]());
        finish ateach(p1 in atoms) {
            val localAtoms = atoms(p1);
            finish for (i in 0..(localAtoms.size-1)) {
                val atom = localAtoms(i);
                val charge = atom.charge;
                val offsetCentre = atom.centre + offset;
                val boxIndex = Fmm3d.getLowestLevelBoxIndex(offsetCentre, lowestLevelDim, size);
                async at(boxAtomsTemp.dist(boxIndex)) {
                    val remoteAtom = new PointCharge(offsetCentre, charge);
                    atomic boxAtomsTemp(boxIndex).add(remoteAtom);
                }
            }
        }

        finish ateach(boxIndex in lowestLevelBoxes) {
            val boxAtoms = boxAtomsTemp(boxIndex);
            if (boxAtoms.size() == 0) {
                // post-prune leaf boxes
                // TODO prune intermediate empty boxes as well
                lowestLevelBoxes(boxIndex) = null;
            } else {
                val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
                box.setAtoms(boxAtoms.toArray());
            }
        }
        timer.stop(TIMER_INDEX_TREE);
    }

    protected def upwardPass() {
        multipoleLowestLevel();
        combineMultipoles();
    }        

    /** 
     * For each atom, creates the multipole expansion and adds it to the
     * lowest level box that contains the atom.
     */
    def multipoleLowestLevel() {
        //Console.OUT.println("multipole lowest level");
        timer.start(TIMER_INDEX_MULTIPOLE);
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val numTerms = this.numTerms; // TODO shouldn't be necessary XTENLANG-1913
        val numLevels = this.numLevels;
        val lowestLevelBoxes = boxes(numLevels);
        finish ateach(p1 in Dist.makeUnique(lowestLevelBoxes.dist.places())) {
            for ([x,y,z] in lowestLevelBoxes.dist(here)) {
                val leafBox = lowestLevelBoxes(x,y,z) as FmmLeafBox;
                if (leafBox != null) {
                    val boxCentre = leafBox.getCentre(size);
                    val boxAtoms = leafBox.getAtoms();
                    for (i in 0..(boxAtoms.size-1)) {
                        val atom = boxAtoms(i);
                        val atomLocation = boxCentre.vector(atom.centre);
                        val atomExpansion = MultipoleExpansion.getOlm(atom.charge, atomLocation, numTerms);
                        leafBox.multipoleExp.unsafeAdd(atomExpansion); // only one thread per box, so unsafeAdd is OK
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_MULTIPOLE);
    }

    /** 
     * For each level above the bottom level, combines multipole expansions for <= 8 child
     * boxes into a single multipole expansion for the parent box.
     */
    def combineMultipoles() {
        timer.start(TIMER_INDEX_COMBINE);
        val fmmOperators = this.fmmOperators; // TODO shouldn't be necessary XTENLANG-1913
        val locallyEssentialTree = this.locallyEssentialTree; // TODO shouldn't be necessary XTENLANG-1913
        for (var level: Int = numLevels-1; level >= topLevel; level--) {
            val thisLevel = level;
            //Console.OUT.println("combine level " + level + " => " + (level-1));
            val thisLevelBoxes = boxes(thisLevel);
            val lowerLevelBoxes = boxes(thisLevel+1);
            val halfSideLength = size / Math.pow2(thisLevel+2);
            finish ateach(p1 in Dist.makeUnique(thisLevelBoxes.dist.places())) {
                val wignerA = fmmOperators().wignerA;
                val complexK = fmmOperators().complexK;
                // can fetch multipoles for lower level while calculating this level
                async Fmm3d.prefetchMultipoles(thisLevel+1, locallyEssentialTree(), lowerLevelBoxes);

                val childMultipoles = new Array[MultipoleExpansion](8);
                for ([x,y,z] in thisLevelBoxes.dist(here)) {
                    val parent = thisLevelBoxes(x,y,z);
                    val childRegion = (2*x)..(2*x+1) * (2*y)..(2*y+1) * (2*z)..(2*z+1);
                    // asynchronously fetch all child multipole expansions
                    finish {
                        for ([x2,y2,z2] in childRegion) async {
                            childMultipoles(childRegion.indexOf(x2, y2, z2)) = getMultipoleExpansionLocalCopy(lowerLevelBoxes, x2, y2, z2);
                        }
                    }
                    // and then sequentially sum them together
                    for ([x2,y2,z2] in childRegion) {
                        val dx = ((x2+1)%2)*2-1;
                        val dy = ((y2+1)%2)*2-1;
                        val dz = ((z2+1)%2)*2-1;
                        val childExp = childMultipoles(childRegion.indexOf(x2, y2, z2));
                        if (childExp != null) {
        	                parent.multipoleExp.translateAndAddMultipole(
			                    Vector3d(dx*halfSideLength, dy*halfSideLength, dz*halfSideLength),
			                    complexK(dx,dy,dz), childExp, wignerA((dx+1)/2, (dy+1)/2, (dz+1)/2));
                        }
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_COMBINE);
    }

    protected def downwardPass():Double {
        timer.start(TIMER_INDEX_DOWNWARD);

        val startingLevel = periodic ? topLevel+1 : topLevel;
        val topLevelExp = periodic ? boxes(0)(0,0,0).localExp : null;
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val boxes = this.boxes; // TODO shouldn't be necessary XTENLANG-1913
        val fmmOperators = this.fmmOperators; // TODO shouldn't be necessary XTENLANG-1913
        val locallyEssentialTree = this.locallyEssentialTree; // TODO shouldn't be necessary XTENLANG-1913
        if (!periodic) {
            // top level boxes won't yet have been fetched
            finish ateach(p1 in Dist.makeUnique()) {
                Fmm3d.prefetchMultipoles(startingLevel, locallyEssentialTree(), boxes(startingLevel));
            }
        }
        val topLevelBoxes = boxes(startingLevel);
        val farField = finish(SumReducer()) {
            ateach([x,y,z] in topLevelBoxes) {
                val box = topLevelBoxes(x,y,z);
                offer box.downward(size, topLevelExp, fmmOperators, locallyEssentialTree, boxes);
            }
        };
        timer.stop(TIMER_INDEX_DOWNWARD);
        return farField / 2.0;
    }

    /**
     * Fetch all multipole expansions required for transform to local.
     * This is communication-intensive, so can be overlapped with computation.
     */
    protected static def prefetchMultipoles(thisLevel:Int, myLET:LocallyEssentialTree, thisLevelBoxes:DistArray[FmmBox](3)) {
        val combinedVList = myLET.combinedVList(thisLevel);
        val thisLevelCopies = myLET.multipoleCopies(thisLevel);
        val vListPlaces = new HashMap[Int,ArrayList[Point(3)]](26); // a place may have up to 26 immediate neighbours

        // separate the vList into partial lists stored at each nearby place
        for (p in combinedVList) {
            val boxIndex = combinedVList(p);
            val placeId = thisLevelBoxes.dist(boxIndex).id;
            var vListForPlace : ArrayList[Point(3)] = vListPlaces.getOrElse(placeId, null);
            if (vListForPlace == null) {
                vListForPlace = new ArrayList[Point(3)]();
                vListPlaces.put(placeId, vListForPlace);
            }
            vListForPlace.add(boxIndex);
        }

        // retrieve the partial list for each place and store into my LET
        finish for(placeEntry in vListPlaces.entries()) async {
            val placeId = placeEntry.getKey();
            val vListForPlace = placeEntry.getValue();
            val vListArray = vListForPlace.toArray();
            val multipolesForPlace = (placeId == here.id) ?
                Fmm3d.getMultipolesForBoxList(thisLevelBoxes, vListArray) :
                at(Place.place(placeId)) { Fmm3d.getMultipolesForBoxList(thisLevelBoxes, vListArray)};
            for (i in 0..(vListArray.size-1)) {
                thisLevelCopies(vListArray(i)) = multipolesForPlace(i);
            }
        }
    }

    /**
     * Given a level and a list of box indexes (as Point(3)) stored
     * at a single place, returns an Array, each element of which 
     * is in turn a MultipoleExpansion for the box.
     */
    private static def getMultipolesForBoxList(thisLevelBoxes : DistArray[FmmBox](3), boxList : Rail[Point(3)]) {
        val multipoleList = new Array[MultipoleExpansion](boxList.size, 
                                                            (i : Int) => 
                                                                getMultipoleExpansionLocalCopy(thisLevelBoxes,
                                                                                     boxList(i)(0), 
                                                                                     boxList(i)(1),
                                                                                     boxList(i)(2))
                                                            );
        return multipoleList;
    }

    /**
     * Fetch all atoms required for direct calculations at each place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchRemoteAtoms() : void {
        timer.start(TIMER_INDEX_PREFETCH);
        val lowestLevelBoxes = boxes(numLevels);
        val locallyEssentialTree = this.locallyEssentialTree; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach(p1 in Dist.makeUnique()) {
            val myLET = locallyEssentialTree();
            val myCombinedUList = myLET.combinedUList;
            val cachedAtoms = myLET.cachedAtoms;

            val uListPlaces = new HashMap[Int,ArrayList[Point(3)]](26); // a place may have up to 26 immediate neighbours
            
            // separate the uList into partial lists stored at each nearby place
            for ([p] in myCombinedUList) {
                val boxIndex = myCombinedUList(p);
                val placeId = lowestLevelBoxes.dist(boxIndex).id;
                var uListForPlace : ArrayList[Point(3)] = uListPlaces.getOrElse(placeId, null);
                if (uListForPlace == null) {
                    uListForPlace = new ArrayList[Point(3)]();
                    uListPlaces.put(placeId, uListForPlace);
                }
                uListForPlace.add(boxIndex);
            }

            // retrieve the partial list for each place and store into my LET
            finish for (placeEntry in uListPlaces.entries()) async {
                val placeId = placeEntry.getKey();
                val uListForPlace = placeEntry.getValue();
                val uListArray = uListForPlace.toArray();
                val atomsForPlace = (placeId == here.id) ?
                    Fmm3d.getAtomsForBoxList(lowestLevelBoxes, uListArray) :
                    at(Place.place(placeId)) { Fmm3d.getAtomsForBoxList(lowestLevelBoxes, uListArray)};
                for (i in 0..(uListArray.size-1)) {
                    myLET.cachedAtoms(uListArray(i)) = atomsForPlace(i);
                }
            }
        }
        timer.stop(TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of box indices as Point(3) stored at a single
     * place, returns an Rail, each element of which is in turn
     * a Rail[PointCharge] containing the atoms for each box.
     */
    private static def getAtomsForBoxList(lowestLevelBoxes : DistArray[FmmBox](3), boxList : Rail[Point(3)]) {
        val atomList = new Rail[Rail[PointCharge]](boxList.size, 
                                                    (i : Int) => 
                                                        getAtomsForBox(lowestLevelBoxes,
                                                                             boxList(i)(0), 
                                                                             boxList(i)(1),
                                                                             boxList(i)(2))
                                                    );
        return atomList;
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

        val lowestLevelBoxes = boxes(numLevels);
        val locallyEssentialTree = this.locallyEssentialTree; // TODO shouldn't be necessary XTENLANG-1913
        val lowestLevelDim = this.lowestLevelDim; // TODO shouldn't be necessary XTENLANG-1913
        val size = this.size; // TODO shouldn't be necessary XTENLANG-1913
        val directEnergy = finish (SumReducer()) {
            ateach(p1 in Dist.makeUnique()) {
                val myLET = locallyEssentialTree();
                val cachedAtoms = myLET.cachedAtoms;
                var thisPlaceEnergy : Double = 0.0;
                for ([x1,y1,z1] in lowestLevelBoxes.dist(here)) {
                    val box1 = lowestLevelBoxes(x1,y1,z1) as FmmLeafBox;
                    if (box1 != null) {
                        val box1Atoms = box1.getAtoms();
                        for (atomIndex1 in 0..(box1Atoms.size-1)) {
                            // direct calculation with all atoms in same box
                            val atom1 = box1Atoms(atomIndex1);
                            for (sameBoxAtomIndex in 0..(atomIndex1-1)) {
                                val sameBoxAtom = box1Atoms(sameBoxAtomIndex);
                                val pairEnergy = atom1.charge * sameBoxAtom.charge / atom1.centre.distance(sameBoxAtom.centre);
                                thisPlaceEnergy += pairEnergy;
                            }
                        }

                        // direct calculation with all atoms in non-well-separated boxes
                        val uList = box1.getUList();
                        for (p in 0..(uList.size-1)) {
                            val boxIndex2 = uList(p);
                            // TODO - should be able to detect Point rank and inline
                            val x2 = boxIndex2(0);
                            val y2 = boxIndex2(1);
                            val z2 = boxIndex2(2);
                            val box2Atoms = cachedAtoms(x2, y2, z2);
                            if (box2Atoms != null) {
                                for (otherBoxAtomIndex in 0..(box2Atoms.size-1)) {
                                    val atom2 = box2Atoms(otherBoxAtomIndex);
                                    for (atomIndex1 in 0..(box1Atoms.size-1)) {
                                        val atom1 = box1Atoms(atomIndex1);
                                        thisPlaceEnergy += atom1.charge * atom2.charge / atom1.centre.distance(atom2.centre);
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

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    private def constructTree(numLevels : Int, topLevel : Int, numTerms : Int, ws : Int, periodic : boolean) 
      : Rail[DistArray[FmmBox](3)] {
        val boxArray = new Array[DistArray[FmmBox](3)](numLevels+1);
        for (thisLevel in topLevel..numLevels) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelDist = MortonDist.makeMorton(0..(levelDim-1) * 0..(levelDim-1) * 0..(levelDim-1));
            if (periodic) {
                boxArray(thisLevel) = DistArray.make[FmmBox](new PeriodicDist(thisLevelDist));
            } else {
                boxArray(thisLevel) = DistArray.make[FmmBox](thisLevelDist);
            }
            //Console.OUT.println("level " + thisLevel + " dist: " + boxArray(thisLevel).dist);
        }

        for (thisLevel in topLevel..(numLevels-1)) {
            val levelDim = Math.pow2(thisLevel) as Int;
            val thisLevelBoxes = boxArray(thisLevel);
            finish ateach([x,y,z] in thisLevelBoxes) {
                val box = new FmmBox(thisLevel, x, y, z, numTerms, Fmm3d.getParentForChild(boxArray, thisLevel, topLevel, x,y,z));
                if (periodic) {
                    if (thisLevel > topLevel) {
                        box.createVListPeriodic(ws);
                    }
                } else {
                    box.createVList(ws);
                }
                thisLevelBoxes(x,y,z) = box;
            }
        }

        val lowestLevelBoxes = boxArray(numLevels);
        finish ateach([x,y,z] in lowestLevelBoxes) {
            val box = new FmmLeafBox(numLevels, x, y, z, numTerms, Fmm3d.getParentForChild(boxArray, numLevels, topLevel, x,y,z));
            if (periodic) {
                box.createUListPeriodic(ws);
                box.createVListPeriodic(ws);
            } else {
                box.createUList(ws);
                box.createVList(ws);
            }
            lowestLevelBoxes(x,y,z) = box;
        }

        return boxArray;
    }

    /**
     * Creates the locally essential tree at the current place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private static def createLocallyEssentialTree(numLevels : Int, topLevel : Int, boxes : Rail[DistArray[FmmBox](3)], periodic : Boolean) : LocallyEssentialTree {
        val lowestLevelBoxes = boxes(numLevels);
        val uMin = new Array[Int](3, Int.MAX_VALUE);
        val uMax = new Array[Int](3, Int.MIN_VALUE);
        val combinedUSet = new HashSet[Point(3)]();
        for ([x,y,z] in lowestLevelBoxes.dist(here)) {
            val box1 = lowestLevelBoxes(x,y,z) as FmmLeafBox;
            if (box1 != null) {
                val uList = box1.getUList();
                for ([p] in uList) {
                    val boxIndex2 = uList(p);
                    if (combinedUSet.add(boxIndex2)) {
                        for (i in 0..2) {
                            uMin(i) = Math.min(uMin(i), boxIndex2(i));
                            uMax(i) = Math.max(uMax(i), boxIndex2(i));        
                        }
                    }
                }
            }
        }
        val combinedUList = new Array[Point(3)](combinedUSet.size());
        var j : Int = 0;
        for (boxIndex in combinedUSet) {
            combinedUList(j++) = boxIndex;
        }

        val combinedVList = new Rail[Rail[Point(3)]](numLevels+1);
        val vListMin = new Rail[Rail[Int]](numLevels+1);
        val vListMax = new Rail[Rail[Int]](numLevels+1);
        for (thisLevel in topLevel..numLevels) {
            if (!(periodic && thisLevel == topLevel)) {
                val vMin = new Array[Int](3, Int.MAX_VALUE);
                val vMax = new Array[Int](3, Int.MIN_VALUE);
                val combinedVSet = new HashSet[Point(3)]();
                val thisLevelBoxes = boxes(thisLevel);
                for ([x,y,z] in thisLevelBoxes.dist(here)) {
                    val box1 = thisLevelBoxes(x,y,z);
                    if (box1 != null) {
                        val vList = box1.getVList();
                        if (vList != null) {
                            for ([p] in vList) {
                                val boxIndex2 = vList(p);
                                if (combinedVSet.add(boxIndex2)) {
                                    for (i in 0..2) {
                                        vMin(i) = Math.min(vMin(i), boxIndex2(i));
                                        vMax(i) = Math.max(vMax(i), boxIndex2(i));        
                                    }
                                }
                            }
                        }
                    }
                }
                //Console.OUT.println("done " + combinedVSet.size());
                vListMin(thisLevel) = vMin;
                vListMax(thisLevel) = vMax;
                val thisLevelVList = new Array[Point(3)](combinedVSet.size());
                var i : Int = 0;
                for (boxIndex in combinedVSet) {
                    thisLevelVList(i++) = boxIndex;
                }
                combinedVList(thisLevel) = thisLevelVList;
            }
        }
        return new LocallyEssentialTree(combinedUList,
                                        combinedVList,
                                        uMin,
                                        uMax,
                                        vListMin,
                                        vListMax);
    }

    protected static def getLowestLevelBoxIndex(offsetCentre : Point3d, lowestLevelDim : Int, size : Double) : Point(3) {
        return  Point.make((offsetCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as Int, (offsetCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as Int);
    }

    protected static def getParentForChild(boxes : Rail[DistArray[FmmBox](3)], level : Int, topLevel : Int, x : Int, y : Int, z : Int) : GlobalRef[FmmBox] {
        if (level == topLevel)
            return GlobalRef[FmmBox](null);
        return (at(boxes(level-1).dist(x/2, y/2, z/2)) {GlobalRef[FmmBox](boxes(level-1)(x/2, y/2, z/2))});
    }

    private static def getAtomsForBox(lowestLevelBoxes : DistArray[FmmBox](3), x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox;
        if (box != null) {
            return box.getAtoms();
        } else {
            return null;
        }
    }

    /**
     * @return a local copy at the current place of a box's multipole expansion
     */
    private static def getMultipoleExpansionLocalCopy(thisLevelBoxes : DistArray[FmmBox](3), x : Int, y : Int, z : Int) : MultipoleExpansion {
        val boxHome = thisLevelBoxes.dist(x,y,z);
        if (boxHome == here) {
            return thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null;
        } else {
            return at(boxHome) {thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null};
        }
    }
}

