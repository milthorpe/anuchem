/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.HashSet;
import x10.util.Pair;
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
 * The 3D simulation space is divided in an octree of up to <code>dMax</code> levels.
 * </p> 
 * @see White & Head-Gordon (1994). "Derivation and efficient implementation 
 *  of the fast multipole method". J Chem Phys 101 (8)
 * @see Lashuk et al. (2009). "A massively parallel adaptive fast-multipole 
 *  method on heterogeneous architectures". Proceedings of SC 2009
 * @author milthorpe
 */
public class FastMultipoleMethod {
    /** The maximum number of levels in the octree. */
    public val dMax:UShort;

    /** 
     * Return the top level of boxes actually used in the method.
     * This is 0 for the periodic FMM and 2 for the non-periodic FMM.
     */
    protected val topLevel:UShort;

    /** The number of lowest level boxes along one side of the cube. */
    public val lowestLevelDim:Int;

    /**
     * The side length of the cube. 
     * To ensure balanced rounding errors within the multipole and local 
     * calculations, all force/energy calculations are performed within 
     * an origin-centred cube.
     */
    public val size:Double; 

    /** The number of terms to use in the multipole and local expansions. */
    public val numTerms:Int;

    /** The well-separatedness parameter ws. */
    public val ws:Int;

    /**
     * Are boundary conditions periodic?
     */
    val periodic:boolean;


    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL:Int = 0;
    public static val TIMER_INDEX_PREFETCH:Int = 1;
    public static val TIMER_INDEX_UPWARD:Int = 2;
    public static val TIMER_INDEX_DOWNWARD:Int = 3;
    public static val TIMER_INDEX_TREE:Int = 4;
    public static val TIMER_INDEX_PLACEHOLDER:Int = 5;

    val localData:PlaceLocalHandle[LocalData];

    static class LocalData {
        /** All leaf octants held at this place. */
        var leafOctants:ArrayList[LeafOctant];

        /** All top-level octants held at this place. */
        var topLevelOctants:ArrayList[Octant];

        /** 
         * The locally essential tree at this place. 
         * @see Lashuk et al. (2009).
         */
        var locallyEssentialTree:LET;

        /** The operator arrays for transformations and translations. */
        val fmmOperators:FmmOperators;

        /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
        public val timer:Timer;

        public def this(numTerms:Int, ws:Int) {
            fmmOperators = new FmmOperators(numTerms, ws);
            timer = new Timer(6);
            // TODO construct LET
        }

    }

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level box
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     */
    public def this(density:Double, 
                    dMax:Int,
                    numTerms:Int,
                    ws:Int,
                    size:Double,  
                    numAtoms:Int) {
        // topLevel in regular FMM is 2 (boxes higher than this cannot be well-spaced)
        this(density, dMax, numTerms, ws, size, numAtoms, 2, false);
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
    protected def this(density:Double, 
                    dMax:Int,
                    numTerms:Int,
                    ws:Int,
                    size:Double,  
                    numAtoms:Int,
                    topLevel:Int,
                    periodic:boolean) {
        this.topLevel = topLevel as UShort;
        this.periodic = periodic;
        this.dMax = dMax as UShort;

        val lowestLevelDim = Math.pow2(dMax);
        this.lowestLevelDim = lowestLevelDim;

        this.numTerms = numTerms;
        this.ws = ws;

        this.size = size;

        this.localData = PlaceLocalHandle.make[LocalData](Dist.makeUnique(), () => new LocalData(numTerms, ws));
    }
    
    public def calculateEnergy():Double {
        val totalEnergy = finish (SumReducer()) {
            ateach(p1 in Dist.makeUnique()) {
                offer calculateEnergyLocal();
            }
        };

        return totalEnergy;
    }

    public def calculateEnergyLocal():Double {
        val timer = localData().timer;
        timer.start(TIMER_INDEX_TOTAL);
        finish {
            async {
//                prefetchRemoteAtoms();
            }
            upwardPass();
            Team.WORLD.barrier(here.id);
        }
        val localEnergy = 0.5 * downwardPass();
        timer.stop(TIMER_INDEX_TOTAL);
        Team.WORLD.allreduce[Long](here.id, timer.total, 0, timer.total, 0, timer.total.size, Team.MAX);

        return localEnergy;

    }

    public def printForces() {
        for(p1 in Place.places()) at(p1) {
            val leafOctants = localData().leafOctants;
            var maxForceError:Double = 0.0;
            for (leafOctant in leafOctants) {
                val atoms = leafOctant.getAtoms();
                for (i in 0..(atoms.size-1)) {
                    val atom = atoms(i);
                    var directForce:Vector3d=Vector3d.NULL;
                    for (j in 0..(atoms.size-1)) {
                        if (i!=j) {
                            val atomJ = atoms(j);
                            val rVec = atomJ.centre - atom.centre;
                            val r2 = rVec.lengthSquared();
                            val r = Math.sqrt(r2);
                            val pairForce = (atom.charge * atomJ.charge / r2 / r) * rVec;
                            directForce += pairForce; 
                        }
                    }
                    for (leafOctant2 in leafOctants) {
                        if (leafOctant2 != leafOctant) {
                            val atoms2 = leafOctant2.getAtoms();
                            for (j in 0..(atoms2.size-1)) {
                                val atomJ = atoms2(j);
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
            Console.OUT.println("max force error = " + maxForceError);
        }
    }

    protected def upwardPass() {
        val local = localData();
        local.timer.start(TIMER_INDEX_UPWARD);
        finish for(octant in local.topLevelOctants) async {
            octant.upward(size, local.fmmOperators, local.locallyEssentialTree, periodic);
        }
        local.timer.stop(TIMER_INDEX_UPWARD);
    }

    protected def downwardPass():Double {
        val local = localData();
        local.timer.start(TIMER_INDEX_DOWNWARD);

        val topLevelExp = null; // TODO periodic ? (at (boxes(0).dist(0,0,0)) {boxes(0)(0,0,0).localExp}) : null;
        val farField = finish(SumReducer()) {
            for(octant in local.topLevelOctants) async {
                offer octant.downward(size, topLevelExp, local.fmmOperators, local.locallyEssentialTree, dMax, periodic);
            }
        };
        local.timer.stop(TIMER_INDEX_DOWNWARD);
        return farField;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    public def initialAssignment(atoms:DistArray[Rail[MMAtom]](1)) {
        finish ateach(p1 in atoms) {
            localData().timer.start(TIMER_INDEX_TREE);
            val localAtoms = atoms(p1);
            assignAtomsToBoxesLocal(localAtoms);
            localData().timer.stop(TIMER_INDEX_TREE);
        }
    }

    /** 
     * As atoms move within the simulation box, their owning place
     * will change.  Send any atoms held at this place to their new
     * owning place.
     */
    public def reassignAtoms(step:Int) {
        val localAtoms = new ArrayList[MMAtom]();
        val leafOctants = localData().leafOctants;
        var numAtoms:Int = 0;
        for (leafOctant in leafOctants) {
            val atoms = leafOctant.getAtoms();
            for (i in atoms) {
                val atom = atoms(i);
                if (atom != null) {
                    numAtoms++;
                    localAtoms.add(atom);
                }
            }
        }

        assignAtomsToBoxesLocal(localAtoms.toArray());

        Team.WORLD.barrier(here.id);

    }

    public def assignAtomsToBoxesLocal(localAtoms:Rail[MMAtom]) {
        val octants = new HashMap[OctantId,LeafOctant]();
        for (i in localAtoms) {
            val atom = localAtoms(i);
            val octantId = getLeafOctantId(atom.centre);
            var octant:LeafOctant = octants.getOrElse(octantId, null);
            if (octant == null) {
                octant = new LeafOctant(octantId, numTerms);
                octants.put(octantId, octant);
            }
            octant.atomList.add(atom);
        }

        val octantEntries = octants.entries();
        val leafOctantList = new ArrayList[LeafOctant]();
        val octantList = new ArrayList[Octant]();
        for (octantEntry in octantEntries) {
            val octant = octantEntry.getValue();
            leafOctantList.add(octant);
            octantList.add(octant);
        }
        leafOctantList.sort();
        val firstId = leafOctantList.getFirst().id;
        val lastId = leafOctantList.getLast().id;
        Console.OUT.println("at " + here + " first leaf = " + firstId + " last = " + lastId);
        localData().leafOctants = leafOctantList;

        // set atom lists for leaves
        //Console.OUT.println("octants at this place:");
        for(octant in leafOctantList) {
            //Console.OUT.println(octant);
            octant.setAtoms(octant.atomList.toArray());
            octant.atomList = new ArrayList[MMAtom](); // clear for next iteration
        }

        // TODO sort and distribute octants

        // TODO coarsening and balancing

        // create shared octants in higher levels
        var level:UShort=dMax-1;
        while (level >= topLevel) {
            val sharedOctants = new HashMap[OctantId,SharedOctant]();
            for (octant in octantList) {
                val parentId = octant.id.getParentId(dMax);
                var parentOctant:SharedOctant = sharedOctants.getOrElse(parentId, null) as SharedOctant;
                if (parentOctant == null) {
                    parentOctant = new SharedOctant(parentId, numTerms);
                    sharedOctants.put(parentId, parentOctant);
                }
                octant.parent = parentOctant;
                parentOctant.children(parentId.getChildIndex(dMax, octant.id)) = octant;
            }

            octantList.clear();
            val sharedOctantEntries = sharedOctants.entries();
            val sharedOctantList = new ArrayList[Octant]();
            for (sharedOctantEntry in sharedOctantEntries) {
                octantList.add(sharedOctantEntry.getValue());
            }
            octantList.sort(); // TODO remove - only useful for logging
            //Console.OUT.println("shared octants for level " + level + ":");
            //for (sharedOctant in octantList) {
            //    Console.OUT.println(sharedOctant);
            //}
            level--;
        }
        localData().topLevelOctants = octantList;

        createLET();
    }

    /**
     * Creates the locally essential tree at the current place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private def createLET() {
        val uMin = new Array[Int](3, Int.MAX_VALUE);
        val uMax = new Array[Int](3, Int.MIN_VALUE);
        val combinedUSet = new HashSet[OctantId]();
        for(octant in localData().leafOctants) {
            octant.createUList(ws);
            val uList = octant.getUList();
            for ([p] in uList) {
                val adjacentOctantId = uList(p);
                combinedUSet.add(adjacentOctantId);
            }
        }

        Console.OUT.println("at " + here + " combined U-list:");
        val combinedUList = new Array[OctantId](combinedUSet.size());
        var j : Int = 0;
        for (octantId in combinedUSet) {
            combinedUList(j++) = octantId;
            Console.OUT.println(octantId);
        }

        val combinedVSet = new HashSet[OctantId]();
        for(octant in localData().topLevelOctants) {
            octant.addToCombinedVSet(combinedVSet, ws);
        }
        //Console.OUT.println("done " + combinedVSet.size());

        Console.OUT.println("at " + here + " combined V-list:");
        val combinedVList = new Rail[OctantId](combinedVSet.size());
        var i : Int = 0;
        for (octantId in combinedVSet) {
            combinedVList(i++) = octantId;
            Console.OUT.println(octantId);
        }
        localData().locallyEssentialTree = new LET(combinedUList, combinedVList);
    }

    /** @return the octant ID for the leaf octant at the maximum depth in the tree */
    public def getLeafOctantId(atomCentre:Point3d):OctantId {
        return new OctantId((atomCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as UShort, (atomCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as UShort, (atomCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as UShort, dMax as UShort);
    }
}

