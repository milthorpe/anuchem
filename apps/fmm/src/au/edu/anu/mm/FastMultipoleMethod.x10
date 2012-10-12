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
import x10.util.ArrayUtils;
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
    public val dMax:UByte;

    /** 
     * Return the top level of octants actually used in the method.
     * This is 0 for the periodic FMM and 2 for the non-periodic FMM.
     */
    protected val topLevel:UByte;

    /** The number of lowest level octants along one side of the cube. */
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

    public static val localData = new FmmLocalData();

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level octant
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
        // topLevel in regular FMM is 2 (octants higher than this cannot be well-spaced)
        this(density, dMax, numTerms, ws, size, numAtoms, 2, false);
    }

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level octant
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     * @param topLevel the topmost level for which octant in the octree are used
     */
    protected def this(density:Double, 
                    dMax:Int,
                    numTerms:Int,
                    ws:Int,
                    size:Double,  
                    numAtoms:Int,
                    topLevel:Int,
                    periodic:boolean) {
        this.topLevel = topLevel as UByte;
        this.periodic = periodic;
        this.dMax = dMax as UByte;

        val lowestLevelDim = Math.pow2(dMax);
        this.lowestLevelDim = lowestLevelDim;

        this.numTerms = numTerms;
        this.ws = ws;

        this.size = size;
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
        val timer = FastMultipoleMethod.localData.timer;
        timer.start(FmmLocalData.TIMER_INDEX_TOTAL);

        finish {
            async {
                prefetchRemoteAtoms();
            }
            upwardPass();
            Team.WORLD.barrier(here.id);

        }
        val potential = downwardPass();
        timer.stop(FmmLocalData.TIMER_INDEX_TOTAL);
        Console.OUT.printf("at %d: prefetch %.3G upward %.3G downward %.3G P2P %.3G\n", here.id, timer.mean(FmmLocalData.TIMER_INDEX_PREFETCH) / 1e9, timer.mean(FmmLocalData.TIMER_INDEX_UPWARD) / 1e9, timer.mean(FmmLocalData.TIMER_INDEX_DOWNWARD) / 1e9, timer.total(FmmLocalData.TIMER_INDEX_P2P) / 1e9);
        Team.WORLD.allreduce[Long](here.id, timer.total, 0, timer.total, 0, timer.total.size, Team.MAX);

        return 0.5 * potential;

    }

    public def printForces() {
        for(p1 in Place.places()) at(p1) {
            val leafOctants = FastMultipoleMethod.localData.leafOctants;
            var maxForceError:Double = 0.0;
            for (leafOctant in leafOctants) {
                val atoms = leafOctant.atoms;
                for (atomI in atoms) {
                    var directForce:Vector3d=Vector3d.NULL;
                    for (atomJ in atoms) {
                        if (atomI != atomJ) {
                            val rVec = atomJ.centre - atomI.centre;
                            val r2 = rVec.lengthSquared();
                            val r = Math.sqrt(r2);
                            val pairForce = (atomI.charge * atomJ.charge / r2 / r) * rVec;
                            directForce += pairForce; 
                        }
                    }
                    for (leafOctant2 in leafOctants) {
                        if (leafOctant2 != leafOctant) {
                            for (atomJ in leafOctant2.atoms) {
                                val rVec = atomJ.centre - atomI.centre;
                                val r2 = rVec.lengthSquared();
                                val r = Math.sqrt(r2);
                                val pairForce = (atomI.charge * atomJ.charge / r2 / r) * rVec;
                                directForce += pairForce;
                            }
                        }
                    }
                    val forceError = (directForce - atomI.force).magnitude() / directForce.magnitude();
                    maxForceError = Math.max(forceError, maxForceError);
                 
                    //Console.OUT.println(atomI.symbol + " force = " + atomI.force + " magnitude " + atomI.force.length() + " forceError = " + forceError);
                }
            }
            Console.OUT.printf("F err: %.2G", maxForceError);
        }
    }

    protected def upwardPass() {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_UPWARD);
        finish for (topLevelOctant in local.topLevelOctants) async {
            topLevelOctant.upward();
        }
        local.timer.stop(FmmLocalData.TIMER_INDEX_UPWARD);
    }

    protected def downwardPass() {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_DOWNWARD);

        val topLevelExp:LocalExpansion = null; // TODO periodic ? (at (boxes(0).dist(0,0,0)) {boxes(0)(0,0,0).localExp}) : null;
        val farField = finish(SumReducer()) {
            for (topLevelOctant in local.topLevelOctants) {
                if (topLevelOctant != null && topLevelOctant.id.level == OctantId.TOP_LEVEL) async {
                    offer topLevelOctant.downward(topLevelExp);
                }
            }
        };
        local.timer.stop(FmmLocalData.TIMER_INDEX_DOWNWARD);
        return farField;
    }

    public def initialAssignment(atoms:DistArray[Rail[MMAtom]](1)) {
        finish ateach(p1 in atoms) {
            FastMultipoleMethod.localData.timer.start(FmmLocalData.TIMER_INDEX_TREE);
            this.localData.init(numTerms, dMax as UByte, ws, size);
            val localAtoms = atoms(p1);
            assignAtomsToOctantsLocal(localAtoms);
            FastMultipoleMethod.localData.timer.stop(FmmLocalData.TIMER_INDEX_TREE);
        }
    }

    /** 
     * Estimate cost of uList and vList interactions at each place, and then
     * combine with estimates at other places to get the average cost.
     */
    private def estimateCostLocal() {
        val estimate = new Array[Long](2);
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        val numEstimates = Math.min(leafOctants.size(), 8);
        var uListEstimate:Long = 0L;
        var vListEstimate:Long = 0L;
        if (numEstimates > 0) {
            for (i in 0..(numEstimates-1)) {
                if (leafOctants(i) != null) {
                    uListEstimate += leafOctants(i).estimateUListCost();
                }
            }
            uListEstimate /= numEstimates;
            vListEstimate = leafOctants(0).estimateVListCost();

        }

        //Console.OUT.println("u-List estimate at " + here + " = " + uListEstimate);
        estimate(FmmLocalData.ESTIMATE_P2P) = uListEstimate;
        //Console.OUT.println("v-List estimate at " + here + " = " + vListEstimate);
        estimate(FmmLocalData.ESTIMATE_M2L) = vListEstimate;

        Team.WORLD.allreduce[Long](here.id, estimate, 0, estimate, 0, estimate.size, Team.ADD);

        val cost = FastMultipoleMethod.localData.cost;
        for (i in 0..(estimate.size-1)) {
            cost(i) = estimate(i) / Place.numPlaces();
        }
        if (here == Place.FIRST_PLACE) {
            Console.OUT.println("u-List cost " + cost(FmmLocalData.ESTIMATE_P2P));
            Console.OUT.println("v-List cost " + cost(FmmLocalData.ESTIMATE_M2L));
        }
    }

    public def countOctants() {
        val totalOctants = finish (IntSumReducer()) {
            ateach(p1 in Dist.makeUnique()) {
                var octants:Int = 0;
                for (topLevelOctant in FastMultipoleMethod.localData.topLevelOctants) {
                    octants += topLevelOctant.countOctants();
                }
                offer octants;
            }
        };
        val totalGhostOctants = finish (IntSumReducer()) {
            ateach(p1 in Dist.makeUnique()) {
                var ghostOctants:Int = 0;
                for (topLevelOctant in FastMultipoleMethod.localData.topLevelOctants) {
                    ghostOctants += topLevelOctant.ghostOctants();
                }
                offer ghostOctants;
            }
        };
        Console.OUT.println("total: " + totalOctants + " ghost: " + totalGhostOctants);
    }

    /** 
     * As atoms move within the simulation cube, their owning place
     * will change.  Send any atoms held at this place to their new
     * owning place.
     */
    public def reassignAtoms(step:Int) {
        FastMultipoleMethod.localData.timer.start(FmmLocalData.TIMER_INDEX_TREE);

        val localAtoms = new ArrayList[MMAtom]();
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        var numAtoms:Int = 0;
        for (leafOctant in leafOctants) {
            for (atom in leafOctant.atoms) {
                if (atom != null) {
                    numAtoms++;
                    localAtoms.add(atom);
                }
            }
        }
        assignAtomsToOctantsLocal(localAtoms.toArray());

        Team.WORLD.barrier(here.id);
        FastMultipoleMethod.localData.timer.stop(FmmLocalData.TIMER_INDEX_TREE);
    }

    public def assignAtomsToOctantsLocal(localAtoms:Rail[MMAtom]) {
        val local = FastMultipoleMethod.localData;
        val timer = local.timer;
        timer.start(FmmLocalData.TIMER_INDEX_ASSIGN);
        val octants = new HashMap[UInt,Octant]();
        for (i in localAtoms) {
            val atom = localAtoms(i);
            val octantId = getLeafOctantId(atom.centre);
            val mortonId = octantId.getMortonId();
            var octant:LeafOctant = octants.getOrElse(mortonId, null) as LeafOctant;
            if (octant == null) {
                octant = new LeafOctant(octantId, numTerms);
                octants.put(mortonId, octant);

            }
            octant.atoms.add(atom);
        }

        val octantEntries = octants.entries();
        val leafOctantList = new ArrayList[LeafOctant]();
        for (octantEntry in octantEntries) {
            val octant = octantEntry.getValue();
            leafOctantList.add(octant as LeafOctant);
        }
        leafOctantList.sort();

        val numTerms = this.numTerms; // TODO shouldn't be necessary XTENLANG-1913
        val octantLoads = local.octantLoads;
        octantLoads.clear();
        if (leafOctantList.size() > 0) {
            for(octant in leafOctantList) {
                val mortonId = octant.id.getLeafMortonId() as Int;
                octantLoads(mortonId) = octant.atoms.size();
                //Console.OUT.println("at " + here + " octant " + mortonId + "(" + octant.id + ") has " + octantLoads(mortonId));
            }
        }

        local.octants = octants;
        local.leafOctants = leafOctantList;

        timer.stop(FmmLocalData.TIMER_INDEX_ASSIGN);

        estimateCostLocal(); // TODO shouldn't include in tree construction time

        // sort and redistribute octants
        timer.start(FmmLocalData.TIMER_INDEX_REDUCE);
        allReduceOctantLoads(octantLoads);
        val numAtoms = octantLoads.reduce[Int]((a:Int,b:Int)=>a+b, 0);
        val uListCost = local.cost(FmmLocalData.ESTIMATE_P2P);
        val q = numAtoms / octantLoads.size; // average particles per lowest-level octant
        val vListCost = local.cost(FmmLocalData.ESTIMATE_M2L);

        val costs = new Array[Long](octantLoads.size);
        for(i in 0..(octantLoads.size-1)) {
            val mortonId = i as UInt;
            //Console.OUT.println(OctantId.getFromMortonId(mortonId) + " uList " + (octantLoads(i) * LeafOctant.estimateUListSize(mortonId, dMax) * q * uListCost) + " vList " + Octant.estimateVListSize(mortonId, dMax) * vListCost);
            // octantLoad = uListInteractions * costP2P + vListInteractions * costM2L 
            costs(i) = 
                    octantLoads(i) * LeafOctant.estimateUListSize(mortonId, dMax) * q * uListCost 
                  + Octant.estimateVListSize(mortonId, dMax) * vListCost;
        }
        timer.stop(FmmLocalData.TIMER_INDEX_REDUCE);
        timer.start(FmmLocalData.TIMER_INDEX_REDIST);
        val total = costs.reduce[Long]((a:Long,b:Long)=>a+b, 0L);
        val placeShare = total / Place.MAX_PLACES;
        val leftOver = total % Place.MAX_PLACES;

        val firstLeafOctant = local.firstLeafOctant;
        firstLeafOctant(0) = 0U;
        var i:Int=0;
        for (p in 1..Place.MAX_PLACES) {
            val share = p <= leftOver ? placeShare+1 : placeShare;
            var load:Long = 0;
            while(load < share && i < costs.size) {
                load += costs(i++);
                if (load >= share) {
                    break;
                }
            }
            firstLeafOctant(p) = i as UInt;
/*
            if (here == Place.FIRST_PLACE) {
                Console.OUT.println("place " + (p-1) + " first " + firstLeafOctant(p-1) + " last " + (firstLeafOctant(p)-1) + " load " + load);
            }
*/
        }

        finish for (var destPlace:Place=here.next(); destPlace!=here && leafOctantList.size() > 0; destPlace=destPlace.next()) {
            val p = destPlace.id;
            val placeStart = firstLeafOctant(p);
            val placeEnd = firstLeafOctant(p+1);
            //Console.OUT.println("at " + here + " sending octants between " + placeStart + " and " + (placeEnd-1) + " to " + destPlace);

            // redistribute octant atoms
            val atomsToSend = new ArrayList[Pair[OctantId,ArrayList[MMAtom]]]();
            atomic {
                var j:Int = leafOctantList.size()-1;
                while(j >= 0 && leafOctantList(j).id.getLeafMortonId() >= placeEnd) j--;
                //Console.OUT.println("at " + here + " leafOctantList.size() = " + leafOctantList.size() + " j =  " + j);
                while(j >= 0 && leafOctantList(j).id.getLeafMortonId() >= placeStart) {
                    val octant = leafOctantList.removeAt(j);
                    atomsToSend.addBefore(0, new Pair[OctantId,ArrayList[MMAtom]](octant.id, octant.atoms));
                    j--;
                }
            }
            //Console.OUT.println("at " + here + " sending " + atomsToSend.size() + " octants to place " + destPlace);
            if (atomsToSend.size() > 0) {
               at(destPlace) async {
                    val localDest = FastMultipoleMethod.localData;
                    val leafOctantsHere = localDest.leafOctants;
                    atomic {
                        var k:Int = 0;
                        for(octantPair in atomsToSend) {
                            val octantId = octantPair.first;
                            val atoms = octantPair.second;
                            //Console.OUT.println("at " + here + " adding " + atoms.size() + " to octant " + octantId + " " + octantId.getLeafMortonId() + " k = " + k + " leafOctantsHere.size() = " + leafOctantsHere.size());
                            var targetOctant:LeafOctant = null;
                            while(k < leafOctantsHere.size()) {
                                val compareOctant = leafOctantsHere(k);
                                val diff = compareOctant.id.compareTo(octantId);
                                if (diff == 0) {
                                    targetOctant = compareOctant;
                                    //Console.OUT.println("at " + here + " found octant to expand " + targetOctant.id);
                                    break;
                                } else if (diff > 0) {
                                    // found insertion point
                                    //Console.OUT.println("inserting " + octantId + " " + octantId.getLeafMortonId() + " before " + compareOctant.id + " " + compareOctant.id.getLeafMortonId());
                                    break;
                                }
                                k++;
                            }
                            if (targetOctant == null) {
                                // target octant has not yet been created at this place
                                //Console.OUT.println("at " + here + " creating octant " + octantId.getLeafMortonId());
                                targetOctant = new LeafOctant(octantId, numTerms);
                                leafOctantsHere.addBefore(k, targetOctant);
                                localDest.octants.put(octantId.getMortonId(), targetOctant);
                            }
                            for (atom in atoms) {
                                targetOctant.atoms.add(atom);
                            }
                        }
                    }
                }
            }
        }
        Team.WORLD.barrier(here.id);
        timer.stop(FmmLocalData.TIMER_INDEX_REDIST);

        // TODO coarsening and balancing
        timer.start(FmmLocalData.TIMER_INDEX_PARENTS);
        val octantList = new ArrayList[Octant](leafOctantList.size());
        var prevOctantId:OctantId = OctantId(0UY,0UY,0UY,dMax);
        for(octant in leafOctantList) {
             //Console.OUT.println("at " + here + ": " + octant.id);
             octantList.add(octant);
             if (octant.id > prevOctantId.next()) {
                 //Console.OUT.println("at " + here + " skipped octants from " + prevOctantId + " " + prevOctantId.getLeafMortonId() + " to " + octant.id + " " + octant.id.getLeafMortonId());
             }
             prevOctantId = octant.id;
        }

        local.topLevelOctants = new ArrayList[Octant](octantList.size());

        // create parent octants in higher levels
        var level:UByte=dMax;
        while (level > OctantId.TOP_LEVEL) {
            val parentOctants = new HashMap[OctantId,ParentOctant]();
            for (octant in octantList) {
                val parentId = octant.id.getParentId();
                val parentAnchorId = parentId.getAnchor(dMax);
                val placeId = local.getPlaceId(parentAnchorId);
                if (placeId == here.id) {
                    //Console.OUT.println("at " + here + " parent for " + octant.id + " = " + parentId);
                    var parentOctant:ParentOctant = parentOctants.getOrElse(parentId, null) as ParentOctant;
                    if (parentOctant == null) {
                        parentOctant = new ParentOctant(parentId, numTerms, dMax);
                        createGhostChildren(parentOctant, octantLoads, firstLeafOctant);
                        parentOctants.put(parentId, parentOctant);
                        octants.put(parentId.getMortonId(), parentOctant);
                    }
                    octant.parent = parentOctant;
                    val childIndex = parentId.getChildIndex(dMax, octant.id);
                    //Console.OUT.println("at " + here + " childIndex for " + octant.id + " = " + childIndex);
                    parentOctant.children(childIndex) = octant;
                } else {
                    //Console.OUT.println("at " + here + " parent for " + octant.id + " = " + parentId + " anchor " + parentAnchorId + " held at " + placeId);
                    local.topLevelOctants.add(octant);
                }
            }

            octantList.clear();
            val parentOctantEntries = parentOctants.entries();
            val parentOctantList = new ArrayList[Octant]();
            for (parentOctantEntry in parentOctantEntries) {
                octantList.add(parentOctantEntry.getValue());
            }

            if (octantList.size() > 0) {
                octantList.sort();
                //Console.OUT.println("at " + here + " level " + level + " first " + octantList.getFirst().id + " last " + octantList.getLast().id);
                var lastOctantId:OctantId = octantList.getLast().id.next();
                //Console.OUT.println("at " + here + " lastOctantId " + lastOctantId);
                while (local.getPlaceId(lastOctantId.getAnchor(dMax)) == here.id) {
                    //Console.OUT.println("at " + here + " filling in octant " + lastOctantId);
                    val fillerOctant = new ParentOctant(lastOctantId, numTerms, dMax);
                    octantList.add(fillerOctant);
                    octants.put(lastOctantId.getMortonId(), fillerOctant);
                    lastOctantId = lastOctantId.next();
                }
            }

            level--;
            //Console.OUT.println("parent octants for level " + level + ":");
            //for (parent in octantList) {
            //    Console.OUT.println(parent);
            //}
        }

        // add remaining octants to top level
        //Console.OUT.println("at " + here + " topLevelOctants:");
        for (octant in octantList) {
            local.topLevelOctants.add(octant);
            //Console.OUT.println("at " + here + "==>" + octant);
        }
        timer.stop(FmmLocalData.TIMER_INDEX_PARENTS);

        createLET();
    }

    private def allReduceOctantLoads(myOctantLoads:Array[Int]) {
        //Team.WORLD.allreduce[Int](here.id, octantLoads, 0, octantLoads, 0, maxLeafOctants, Team.ADD);
        Team.WORLD.barrier(here.id);

        if (here == Place.FIRST_PLACE) {
            for (p in 1..(Place.MAX_PLACES-1)) {
                val pLoads = at(Place(p)) FastMultipoleMethod.localData.octantLoads;
                myOctantLoads.map(myOctantLoads, pLoads, (x:Int,y:Int)=>x+y);
            }
            finish for(p in 1..(Place.MAX_PLACES-1)) at(Place(p)) async {
                Array.copy[Int](myOctantLoads, FastMultipoleMethod.localData.octantLoads);
            }

        }

        Team.WORLD.barrier(here.id);
    }

    private def createGhostChildren(parentOctant:ParentOctant, octantLoads:Rail[Int], firstLeafOctant:Rail[UInt]) {
        val id = parentOctant.id;
        val levelDim = (Math.pow2(dMax) / Math.pow2(id.level));
        var i:Int = 0;
        for (x2 in (2*id.x)..(2*id.x+1)) {
            for (y2 in (2*id.y)..(2*id.y+1)) {
                for (z2 in (2*id.z)..(2*id.z+1)) {
                    val childOctantId = OctantId(x2 as UByte, y2 as UByte, z2 as UByte, id.level+1U);
                    val anchor = childOctantId.getAnchor(dMax);
                    var nonEmpty:Boolean = false;
                    val size = Math.pow(8.0, (dMax / childOctantId.level - 1)) as Int;
                    val startDescendant = anchor.getLeafMortonId() as Int;
                    val endDescendant = startDescendant + size - 1;
                    for (j in Math.max(firstLeafOctant(here.id+1) as Int,startDescendant)..endDescendant) {
                        if (octantLoads(j) > 0) {
                            nonEmpty = true;
                            break;
                        }
                    }
                    if (nonEmpty) {
                        // non-empty child octant is not held at this place
                        val placeId = FastMultipoleMethod.localData.getPlaceId(anchor);
                        if (placeId != here.id) {
                            //Console.OUT.println("at " + here + " octant " + id + " creating ghost " + childOctantId + " held at " + placeId);
                            parentOctant.children(i) = new GhostOctant(childOctantId, placeId);
                       }
                    }
                    i++;
                }
            }
        }
    }

    /**
     * Creates the locally essential tree at the current place.  This is
     * later used to overlap remote retrieval of multipole expansion and
     * particle data with other computation.
     */
    private def createLET() {
        val local = FastMultipoleMethod.localData;
        val timer = local.timer;
        timer.start(FmmLocalData.TIMER_INDEX_LET);
        val uMin = new Array[Int](3, Int.MAX_VALUE);
        val uMax = new Array[Int](3, Int.MIN_VALUE);
        val combinedUSet = new HashSet[UInt]();
        for(octant in local.leafOctants) {
            octant.createUList(ws);
            val uList = octant.getUList();
            for ([p] in uList) {
                val adjacentOctantMortonId = uList(p).getMortonId();
                combinedUSet.add(adjacentOctantMortonId);
            }
        }

        //Console.OUT.println("at " + here + " combined U-list:");
        val combinedUList = new Array[UInt](combinedUSet.size());
        var j : Int = 0;
        for (mortonId in combinedUSet) {
            combinedUList(j++) = mortonId;
            //Console.OUT.println(mortonId);
        }
        ArrayUtils.sort(combinedUList);

        val combinedVSet = new HashSet[UInt]();
        for (topLevelOctant in local.topLevelOctants) {
            topLevelOctant.addToCombinedVSet(combinedVSet, ws);
        }
        //Console.OUT.println("done " + combinedVSet.size());

        //Console.OUT.println("at " + here + " combined V-list:");
        val combinedVList = new Rail[UInt](combinedVSet.size());
        var i:Int = 0;
        for (mortonId in combinedVSet) {
            combinedVList(i++) = mortonId;
            //Console.OUT.println(mortonId);
        }
        ArrayUtils.sort(combinedVList);
        local.locallyEssentialTree = new LET(combinedUList, combinedVList);
        timer.stop(FmmLocalData.TIMER_INDEX_LET);
    }

    /**
     * Fetch all atoms required for direct calculations at this place.
     * This is communication-intensive, so can be overlapped with computation.
     */
    def prefetchRemoteAtoms() : void {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_PREFETCH);

        val myLET = local.locallyEssentialTree;
        val firstLeafOctant = local.firstLeafOctant;
        val myCombinedUList = myLET.combinedUList;
        val cachedAtoms = myLET.cachedAtoms;

        val uListPlaces = new HashMap[Int,ArrayList[UInt]](26); // a place may have up to 26 immediate neighbours
        
        // separate the uList into partial lists stored at each nearby place
        for ([p] in myCombinedUList) {
            val octantId = OctantId.getFromMortonId(myCombinedUList(p));
            val placeId = local.getPlaceId(octantId.getAnchor(dMax));
            if (placeId >= 0 && placeId < Place.MAX_PLACES) {
                // leaf octant exists at Place(placeId)
                var uListForPlace : ArrayList[UInt] = uListPlaces.getOrElse(placeId, null);
                if (uListForPlace == null) {
                    uListForPlace = new ArrayList[UInt]();
                    uListPlaces.put(placeId, uListForPlace);
                }
                uListForPlace.add(myCombinedUList(p));
            }
        }

        // retrieve the partial list for each place and store into my LET
        finish for (placeEntry in uListPlaces.entries()) async {
            val placeId = placeEntry.getKey();
            val uListForPlace = placeEntry.getValue();
            val uListArray = uListForPlace.toArray();
            val atomsForPlace = (placeId == here.id) ?
                FastMultipoleMethod.getAtomsForOctantList(uListArray) :
                at(Place.place(placeId)) { FastMultipoleMethod.getAtomsForOctantList(uListArray)};
            for (i in 0..(uListArray.size-1)) {
                myLET.setAtomDataForOctant(uListArray(i), atomsForPlace(i));
            }
        }
        local.timer.stop(FmmLocalData.TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of octant indexes stored at a single
     * place, returns an Rail, each element of which is in turn
     * a Rail[PointCharge] containing the atoms for each octant.
     */
    private static def getAtomsForOctantList(octantList:Rail[UInt]) {
        val local = FastMultipoleMethod.localData;
        val atomList = new Rail[Rail[Double]](octantList.size);
        for (i in 0..(octantList.size-1)) {
            val octant = local.getOctant(octantList(i)) as LeafOctant;
            if (octant != null) {
                atomList(i) = octant.getAtomData();
            }
        }
        return atomList;
    }

    /** @return the octant ID for the leaf octant at the maximum depth in the tree */
    public def getLeafOctantId(atomCentre:Point3d):OctantId {
        return new OctantId((atomCentre.i / size * lowestLevelDim + lowestLevelDim / 2) as UByte, (atomCentre.j / size * lowestLevelDim + lowestLevelDim / 2) as UByte, (atomCentre.k / size * lowestLevelDim + lowestLevelDim / 2) as UByte, dMax as UByte);
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    static struct IntSumReducer implements Reducible[Int] {
        public def zero() = 0;
        public operator this(a:Int, b:Int) = (a + b);
    }


}

