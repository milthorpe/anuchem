/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.HashSet;
import x10.util.Team;
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
     * The side length of the cube. 
     * To ensure balanced rounding errors within the multipole and local 
     * calculations, all force/energy calculations are performed within 
     * an origin-centred cube.
     */
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
    public static val TIMER_INDEX_UPWARD : Int = 2;
    public static val TIMER_INDEX_DOWNWARD : Int = 3;
    public static val TIMER_INDEX_TREE : Int = 4;
    public static val TIMER_INDEX_PLACEHOLDER : Int = 5;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = PlaceLocalHandle.make[Timer](PlaceGroup.WORLD, ()=>new Timer(6));

    /** All boxes in the octree division of space. 
     * Array has numLevels elements, for levels [1..numLevels]
     * where boxes(n) = boxes at level n+1
     * Array dimensions are:
     * 0: x coordinate at that level (range 0..2^level)
     * 1: y coordinate
     * 2: z coordinate
     */
    public val boxes : Rail[DistArray[FmmBox](3)];

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
                    size : Double,  
                    numAtoms : Int) {
        // topLevel in regular FMM is 2 (boxes higher than this cannot be well-spaced)
        this(density, numTerms, ws, size, numAtoms, 2, false);
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

        this.size = size;

        timer().start(TIMER_INDEX_TREE);
        val boxes = constructTree(numLevels, topLevel, numTerms, ws, periodic);
        this.boxes = boxes;
        this.locallyEssentialTree = PlaceLocalHandle.make[LocallyEssentialTree](Dist.makeUnique(), () => Fmm3d.createLocallyEssentialTree(numLevels, topLevel, boxes, periodic));

        this.fmmOperators = PlaceLocalHandle.make[FmmOperators](Dist.makeUnique(), () => new FmmOperators(numTerms, ws));

        timer().stop(TIMER_INDEX_TREE);
    }
    
    public def calculateEnergy() : Double {
        timer().start(TIMER_INDEX_TOTAL);
        val totalEnergy = finish (SumReducer()) {
            ateach(p1 in Dist.makeUnique()) {
                finish {
                    async {
                        prefetchRemoteAtoms();
                    }
                    upwardPass();
                    Team.WORLD.barrier(here.id);
                }
                val localEnergy = 0.5 * downwardPass();
                offer localEnergy;
            }
        };

        timer().stop(TIMER_INDEX_TOTAL);

        // reduce timer totals
        finish ateach(p1 in Dist.makeUnique()) {
            Team.WORLD.allreduce[Long](here.id, timer().total, 0, timer().total, 0, timer().total.size, Team.MAX);
        }

        return totalEnergy;
    }

    public def printForces() {
        val lowestLevelBoxes = boxes(numLevels);
        for(p1 in Place.places()) at(p1) {
            var maxForceError:Double = 0.0;
            for ([x,y,z] in lowestLevelBoxes.dist(here)) {
                val leafBox = lowestLevelBoxes(x,y,z) as FmmLeafBox;
                if (leafBox != null) {
                    val boxAtoms = leafBox.getAtoms();
                    for (i in 0..(boxAtoms.size-1)) {
                        val atom = boxAtoms(i);
                        var directForce:Vector3d=Vector3d.NULL;
                        for (j in 0..(boxAtoms.size-1)) {
                            if (i!=j) {
                                val atomJ = boxAtoms(j);
                                val rVec = atomJ.centre - atom.centre;
                                val r2 = rVec.lengthSquared();
                                val r = Math.sqrt(r2);
                                val pairForce = (atom.charge * atomJ.charge / r2 / r) * rVec;
                                directForce += pairForce; 
                            }
                        }
                        for ([x2,y2,z2] in lowestLevelBoxes.dist(here)) {
                            val box2 = lowestLevelBoxes(x2,y2,z2) as FmmLeafBox;
                            if ((x != x2 || y != y2 || z != z2) && box2 != null) {
                                val box2Atoms = box2.getAtoms();
                                for (j in 0..(box2Atoms.size-1)) {
                                    val atomJ = box2Atoms(j);
                                    val rVec = atomJ.centre - atom.centre;
                                    val r2 = rVec.lengthSquared();
                                    val r = Math.sqrt(r2);
                                    val pairForce = (atom.charge * atomJ.charge / r2 / r) * rVec;
                                    directForce += pairForce;
                                }
                            }
                        }
                        val forceError = (directForce - atom.force).magnitude() / directForce.magnitude();
                        maxForceError = Math.max(forceError, maxForceError);
                     
                        //Console.OUT.println(atom.symbol + " force = " + atom.force + " magnitude " + atom.force.length() + " forceError = " + forceError);
                    }
                }
            }
            Console.OUT.println("max force error = " + maxForceError);
        }
    }


    public def assignAtomsToBoxes(atoms:DistArray[Rail[MMAtom]](1)) {
        timer().start(TIMER_INDEX_TREE);
        finish ateach(p1 in atoms) {
            val localAtoms = atoms(p1);
            assignAtomsToBoxesLocal(localAtoms);
            Team.WORLD.barrier(here.id);
            pruneTreeLocal();
        }
        timer().stop(TIMER_INDEX_TREE);
    }

    public def assignAtomsToBoxesLocal(localAtoms:Rail[MMAtom]) {
        val lowestLevelBoxes = boxes(numLevels);
        finish for (i in 0..(localAtoms.size-1)) {
            val atom = localAtoms(i);
            val boxIndex = Fmm3d.getLowestLevelBoxIndex(atom.centre, lowestLevelDim, size);
            val atomPlace = lowestLevelBoxes.dist(boxIndex);
            if (atomPlace == here) {
                val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
                val copyAtom = new MMAtom(atom);
                atomic box.atomList.add(copyAtom);
            } else {
                at(atomPlace) async {
                    val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
                    atomic box.atomList.add(atom);
                }
            }
        }

    }

    public def pruneTreeLocal() {
        val lowestLevelBoxes = boxes(numLevels);
        for(boxIndex in lowestLevelBoxes.dist(here)) {
            val box = lowestLevelBoxes(boxIndex) as FmmLeafBox;
            if (box.atomList.size() == 0) {
                // post-prune leaf boxes
                // TODO prune intermediate empty boxes as well
                lowestLevelBoxes(boxIndex) = null;
            } else {
                box.setAtoms(box.atomList.toArray());
                box.atomList = new ArrayList[MMAtom](); // clear for next iteration
            }
        }
    }

    /** @return a Rail containing all atoms held at this place */
    public def getAtomsLocal() {
        val guessNumAtoms = (Math.pow(8.0, numLevels) * 50) as Int;
        val atomList = new ArrayList[MMAtom](guessNumAtoms);
        val lowestLevelBoxes = boxes(numLevels);
        for([x,y,z] in lowestLevelBoxes.dist(here)) {
            val box = lowestLevelBoxes(x,y,z) as FmmLeafBox;
            if (box != null) {
                val atoms = box.getAtoms();
                for (i in atoms) {
                    val atom = atoms(i);
                    atomList.add(atom);
                }
            }
        }
        return atomList.toArray();
    }

    protected def upwardPass() {
        timer().start(TIMER_INDEX_UPWARD);

        val startingLevel = periodic ? topLevel+1 : topLevel;
        val topLevelBoxes = boxes(startingLevel);
        finish for([x,y,z] in topLevelBoxes.dist(here)) {
            val box = topLevelBoxes(x,y,z);
            if (box != null) {
                async box.upward(size, fmmOperators, locallyEssentialTree, boxes, periodic);
            }
        }
        timer().stop(TIMER_INDEX_UPWARD);
    }

    protected def downwardPass():Double {
        timer().start(TIMER_INDEX_DOWNWARD);

        val startingLevel = periodic ? topLevel+1 : topLevel;
        val topLevelExp = periodic ? (at (boxes(0).dist(0,0,0)) {boxes(0)(0,0,0).localExp}) : null;
        val topLevelBoxes = boxes(startingLevel);

        val farField = finish(SumReducer()) {
            for([x,y,z] in topLevelBoxes.dist(here)) {
                val box = topLevelBoxes(x,y,z);
                if (box != null) {
                    async {
                        offer box.downward(size, topLevelExp, fmmOperators, locallyEssentialTree, boxes, numLevels, periodic);
                    }
                }
            }
        };
        timer().stop(TIMER_INDEX_DOWNWARD);
        return farField;
    }

    /**
     * Fetch all multipole expansions for the combined V-list == boxes that are
     *  well separated from boxes held at this place.
     * @param level the tree level for which multipole expansions are to be fetched
     * @param myLET the LocallyEssentialTree at this place
     * @param boxes the distributed array of boxes at the given level
     */
    protected static def fetchMultipoles(level:Int, myLET:LocallyEssentialTree, boxes:DistArray[FmmBox](3)) {
        val vListPlaces = new HashMap[Int,ArrayList[Point(3)]](26); // a place may have up to 26 immediate neighbours

        val combinedVList = myLET.combinedVList(level);
        // separate the vList into partial lists stored at each nearby place
        if (combinedVList != null) {
            for (p in combinedVList) {
                val boxIndex = combinedVList(p);
                val placeId = boxes.dist(boxIndex).id;
                var vListForPlace : ArrayList[Point(3)] = vListPlaces.getOrElse(0, null);
                if (vListForPlace == null) {
                    vListForPlace = new ArrayList[Point(3)]();
                    vListPlaces.put(placeId, vListForPlace);
                }
                vListForPlace.add(boxIndex);
            }
        }

        val multipoleCopies = myLET.multipoleCopies(level);

        // retrieve the multipole copies from each place and store into my LET
        finish for(placeEntry in vListPlaces.entries()) async {
            val placeId = placeEntry.getKey();
            val vListForPlace = placeEntry.getValue();
            val vListArray = vListForPlace.toArray();
            val multipolesForPlace = (placeId == here.id) ?
                Fmm3d.getMultipolesForBoxList(boxes, vListArray) :
                at(Place.place(placeId)) { Fmm3d.getMultipolesForBoxList(boxes, vListArray)};
            for (i in 0..(vListArray.size-1)) {
                multipoleCopies(vListArray(i)) = multipolesForPlace(i);
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
     * Fetch all atoms required for direct calculations at this place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchRemoteAtoms() : void {
        timer().start(TIMER_INDEX_PREFETCH);
        val lowestLevelBoxes = boxes(numLevels);
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
        timer().stop(TIMER_INDEX_PREFETCH);
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

    public static def getLowestLevelBoxIndex(atomCentre : Point3d, lowestLevelDim : Int, size : Double) : Point(3) {
        return  Point.make((atomCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as Int, (atomCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as Int, (atomCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as Int);
    }

    protected static def getParentForChild(boxes : Rail[DistArray[FmmBox](3)], level : Int, topLevel : Int, x : Int, y : Int, z : Int) : GlobalRef[FmmBox] {
        if (level == topLevel)
            return GlobalRef[FmmBox](null);
        return (at(boxes(level-1).dist(x/2, y/2, z/2)) {GlobalRef[FmmBox](boxes(level-1)(x/2, y/2, z/2))});
    }

    private static def getAtomsForBox(lowestLevelBoxes : DistArray[FmmBox](3), x : Int, y : Int, z : Int) {
        val box = lowestLevelBoxes(x, y, z) as FmmLeafBox;
        if (box != null) {
            return box.getAtomCharges();
        } else {
            return null;
        }
    }

    /**
     * @return a local copy at the current place of a box's multipole expansion
     */
    public static def getMultipoleExpansionLocalCopy(thisLevelBoxes : DistArray[FmmBox](3), x : Int, y : Int, z : Int) : MultipoleExpansion {
        val boxHome = thisLevelBoxes.dist(x,y,z);
        if (boxHome == here) {
            return thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null;
        } else {
            return at(boxHome) {thisLevelBoxes(x,y,z) != null ? (thisLevelBoxes(x,y,z)).multipoleExp : null};
        }
    }
}

