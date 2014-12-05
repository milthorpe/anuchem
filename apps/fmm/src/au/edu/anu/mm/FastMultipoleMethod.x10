/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012-2013.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.regionarray.Dist;
import x10.util.ArrayList;
import x10.util.RailUtils;
import x10.util.HashMap;
import x10.util.Pair;
import x10.util.Team;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.util.Timer;
import au.edu.anu.chem.mm.ParticleData;

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
    /** Estimated density per lowest-level box */
    public val density:Int;
    
    /** The maximum number of levels in the octree. */
    public val dMax:UByte;

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

    /** Verbose logging */
    val verbose:boolean;

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
                    size:Double) {
        this(density, dMax, numTerms, ws, size, false);
    }

    /**
     * Initialises a fast multipole method electrostatics calculation
     * for the given system of atoms.
     * @param density mean number of particles per lowest level octant
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     * @param size length of a side of the simulation cube
     * @param atoms the atoms for which to calculate electrostatics
     */
    protected def this(density:Double, 
                    dMax:Int,
                    numTerms:Int,
                    ws:Int,
                    size:Double,
                    verbose:Boolean
                    ) {
        this.density = density as Int;
        this.dMax = dMax as UByte;

        val lowestLevelDim = Math.pow2(dMax);
        this.lowestLevelDim = lowestLevelDim;

        this.numTerms = numTerms;
        this.ws = ws;

        this.size = size;
        this.verbose = verbose;
    }
    
    public def calculateEnergy():Double {
        val totalEnergy = finish (Reducible.SumReducer[Double]()) {
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
            multipolesToLocal();
        }
        val potential = downwardPass();
        timer.stop(FmmLocalData.TIMER_INDEX_TOTAL);
        if (verbose) {
            Console.OUT.printf("at %d: prefetch %.3G upward %.3G M2L %.3G downward %.3G\n", here.id, timer.mean(FmmLocalData.TIMER_INDEX_PREFETCH) / 1e9, timer.mean(FmmLocalData.TIMER_INDEX_UPWARD) / 1e9, timer.mean(FmmLocalData.TIMER_INDEX_M2L) / 1e9, timer.mean(FmmLocalData.TIMER_INDEX_DOWNWARD) / 1e9);
        }

        return 0.5 * potential;

    }

    public def reduceMaxTimes() {
        finish ateach(p1 in Dist.makeUnique()) {
            val timer = FastMultipoleMethod.localData.timer;
            Team.WORLD.allreduce[Long](timer.total, 0L, timer.total, 0L, timer.total.size, Team.MAX);
        }
    }

    public def printRMSErrors() {
        val particleData = FastMultipoleMethod.localData.particleData;
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        val plm = FmmScratch.getWorkerLocal().plm;
        var mseForce:Double = 0.0;
        var normForce:Double = 0.0;
        var msePot:Double = 0.0;
        var normPot:Double = 0.0;
        var minForce:Double = Double.MAX_VALUE;
        var maxForce:Double = 0.0;
        var minPot:Double = Double.MAX_VALUE;
        var maxPot:Double = 0.0;

        for (leafOctant in leafOctants) {
            val atoms = leafOctant.atoms;
            for (atomI in atoms) {
                val force = particleData.fx(atomI); // because force is overwritten by calculatePotentialAndForces
                var directForce:Vector3d=Vector3d.NULL;
                var directPotential:Double = 0.0;

                val boxCentre = leafOctant.getCentre(size);
                val locationWithinBox = particleData.x(atomI).vector(boxCentre);
                var potential:Double = leafOctant.localExp.calculatePotentialAndForces(atomI, particleData, locationWithinBox, plm);

                val atomICharge = particleData.atomTypes(particleData.atomTypeIndex(atomI)).charge;
                var nonWsPotential:Double = 0.0;
                for (atomJ in atoms) {
                    if (atomI != atomJ) {
                        val rVec = particleData.x(atomJ) - particleData.x(atomI);
                        val invR2 = 1.0 / rVec.lengthSquared();
                        val invR = Math.sqrt(invR2);
                        val e = atomICharge
                              * particleData.atomTypes(particleData.atomTypeIndex(atomJ)).charge
                              * invR;
                        directPotential += e;
                        potential += e;
                        nonWsPotential += e;
                        val pairForce = e * invR2 * rVec;
                        directForce += pairForce; 
                    }
                }

                for (leafOctant2 in leafOctants) {
                    if (leafOctant2 != leafOctant) {
                        for (atomJ in leafOctant2.atoms) {
                            val rVec = particleData.x(atomJ) - particleData.x(atomI);
                            val invR2 = 1.0 / rVec.lengthSquared();
                            val invR = Math.sqrt(invR2);
                            val e = atomICharge
                                  * particleData.atomTypes(particleData.atomTypeIndex(atomJ)).charge
                                  * invR;
                            directPotential += e;
                            val pairForce = e * invR2 * rVec;
                            directForce += pairForce;
                            if (Math.abs(leafOctant2.id.x - leafOctant.id.x as Int) <= 1
                             && Math.abs(leafOctant2.id.y - leafOctant.id.y as Int) <= 1
                             && Math.abs(leafOctant2.id.z - leafOctant.id.z as Int) <= 1) {
                                potential += e;
                                nonWsPotential += e;
                            }
                        }
                    }
                }
                //Console.OUT.println("force = " + force + " directForce = "  + directForce);

                val forceMag = directForce.lengthSquared();
                minForce = Math.min(forceMag, minForce);
                maxForce = Math.max(forceMag, maxForce); 
                mseForce += (force - directForce).lengthSquared();
                normForce += forceMag;

                //Console.OUT.println("potential = " + potential + " directPotential = "  + directPotential);

                val direct = directPotential / atomICharge;
                val fmm = potential / atomICharge;

                minPot = Math.min(direct, minPot);
                maxPot = Math.max(direct, maxPot); 
                msePot += (fmm - direct) * (fmm - direct);
                normPot += direct * direct;
            }
        }

        //Console.OUT.println("minPot = " + minPot + " maxPot = " + maxPot);
        //Console.OUT.println("msePot = " + msePot + " normPot = " + normPot);
        Console.OUT.printf("RMS relative potential err: %.1E\n", Math.sqrt(msePot/normPot));

        //Console.OUT.println("minForce = " + minForce + " maxForce = " + maxForce);
        //Console.OUT.println("mseForce = " + mseForce + " normForce = " + normForce);
        Console.OUT.printf("RMS relative force err: %.1E\n", Math.sqrt(mseForce/normForce));
    }

    protected def upwardPass() {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_UPWARD);
        finish for (topLevelOctant in local.topLevelOctants) async {
            topLevelOctant.upward();
        }
        // force completion of sendMultipoles at all places
        Team.WORLD.barrier();
        local.timer.stop(FmmLocalData.TIMER_INDEX_UPWARD);
    }

    protected def multipolesToLocal() {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_M2L);

        finish for (topLevelOctant in local.topLevelOctants) async {
            topLevelOctant.multipolesToLocal();
        }
        local.timer.stop(FmmLocalData.TIMER_INDEX_M2L);
    }

    protected def downwardPass() {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_DOWNWARD);

        val topLevelExp:LocalExpansion = null; // TODO periodic ? (at (boxes(0).dist(0,0,0)) {boxes(0)(0,0,0).localExp}) : null;
        val farField = finish(Reducible.SumReducer[Double]()) {
            for (topLevelOctant in local.topLevelOctants) {
                if (topLevelOctant != null && topLevelOctant.id.level == OctantId.TOP_LEVEL) async {
                    offer topLevelOctant.downward(topLevelExp);
                }
            }
        };
        local.timer.stop(FmmLocalData.TIMER_INDEX_DOWNWARD);
        return farField;
    }

    public def initialAssignment(numAtoms:Int, particleDataPlh:PlaceLocalHandle[ParticleData]) {
        finish ateach(p1ace in Dist.makeUnique()) {
            FastMultipoleMethod.localData.timer.start(FmmLocalData.TIMER_INDEX_TREE);
            this.localData.init(numTerms, dMax as UByte, ws, size);
            FmmScratch.init( ()=>new FmmScratch(numTerms) );
            estimateCostLocal(numAtoms); // TODO shouldn't include in tree construction time

            val particleData = particleDataPlh();
            FastMultipoleMethod.localData.particleData = particleData;

            // get leaf Morton ID for each atom
            val octantIds = new Rail[UInt](particleData.numAtoms());
            val invSideLength = lowestLevelDim / size;
            val offset = lowestLevelDim / 2.0;
            for (i in 0..(particleData.numAtoms()-1)) {
                octantIds(i) = OctantId.getLeafMortonId(particleData.x(i), invSideLength, offset);
            }
            assignAtomsToOctantsLocal(octantIds);
            FastMultipoleMethod.localData.timer.stop(FmmLocalData.TIMER_INDEX_TREE);
        }
    }

    /** 
     * Estimate cost of uList and vList interactions at each place, and then
     * combine with estimates at other places to get the average cost.
     */
    private def estimateCostLocal(numAtoms:Int) {
        val estimate = new Rail[Long](2);
        val uListEstimate = LeafOctant.estimateUListCost(density);
        val vListEstimate = Octant.estimateVListCost(numTerms);

        estimate(FmmLocalData.ESTIMATE_P2P) = uListEstimate;
        estimate(FmmLocalData.ESTIMATE_M2L) = vListEstimate;

        Team.WORLD.allreduce[Long](estimate, 0L, estimate, 0L, estimate.size, Team.ADD);

        val cost = FastMultipoleMethod.localData.cost;
        for (i in 0..(estimate.size-1)) {
            cost(i) = estimate(i) / Place.numPlaces();
        }
        if (verbose && here == Place.FIRST_PLACE) {
            Console.OUT.println("u-List cost " + cost(FmmLocalData.ESTIMATE_P2P));
            Console.OUT.println("v-List cost " + cost(FmmLocalData.ESTIMATE_M2L));
        }
    }

    public def countOctants() {
        val totalOctants = finish (Reducible.SumReducer[Int]()) {
            ateach(p1 in Dist.makeUnique()) {
                var octants:Int = 0n;
                for (topLevelOctant in FastMultipoleMethod.localData.topLevelOctants) {
                    octants += topLevelOctant.countOctants();
                }
                offer octants;
            }
        };
        val totalGhostOctants = finish (Reducible.SumReducer[Int]()) {
            ateach(p1 in Dist.makeUnique()) {
                var ghostOctants:Int = 0n;
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
    public def reassignAtoms(step:Long) {
        FastMultipoleMethod.localData.timer.start(FmmLocalData.TIMER_INDEX_TREE);
        val particleData = FastMultipoleMethod.localData.particleData;
        val leafOctants = FastMultipoleMethod.localData.leafOctants;
        val octantIds = new Rail[UInt](particleData.numAtoms());
        val invSideLength = lowestLevelDim / size;
        val offset = lowestLevelDim / 2.0;
        for (leafOctant in leafOctants) {
            //Console.OUT.println(here + " leafOctant " + leafOctant.id + " has " + leafOctant.atoms.size() + " atoms");
            for (localAtomIndex in leafOctant.atoms) {
                octantIds(localAtomIndex) = OctantId.getLeafMortonId(particleData.x(localAtomIndex), invSideLength, offset);
                //Console.OUT.println(here + " reassign " + localAtomIndex + " center " + particleData.x(localAtomIndex) + " to " + octantIds(localAtomIndex));
            }
        }
        assignAtomsToOctantsLocal(octantIds);

        Team.WORLD.barrier();
        FastMultipoleMethod.localData.timer.stop(FmmLocalData.TIMER_INDEX_TREE);
    }

    /**
     * Assign all atoms to boxes, redistributing to other places to balance
     * load as necessary.
     * @param localAtoms an unsorted array of atoms with leaf octant Morton IDs
     */
    private def assignAtomsToOctantsLocal(octantIds:Rail[UInt]) {

        val octantAtoms = sortAtoms(octantIds);

        determineLoadBalanceLocal(octantAtoms);

        redistributeOctantsLocal(octantAtoms);

        // TODO coarsening and balancing

        createParentOctantsLocal();

        createLET();
    }

    private def sortAtoms(octantIds:Rail[UInt]) {
        val local = FastMultipoleMethod.localData;
        val timer = local.timer;

        timer.start(FmmLocalData.TIMER_INDEX_SORT);
        val particleData = local.particleData;
        local.octants.clear();

        val numBoxes = Math.pow(8.0, dMax as Int) as Int;

        // histogram sort of atoms
        val octantLoads = local.octantLoads;
        octantLoads.clear();
        var filled:Long = 0;
        for (i in 0..(octantIds.size-1)) {
            val leafMortonId = octantIds(i) as Int;
            //Console.OUT.println(here + " atom " + particleData.globalIndex(i) + " is in " + leafMortonId);
            if (octantLoads(leafMortonId)++ == 0L) filled++;
        }
        //Console.OUT.println(here + " filled " + filled + " octants");
        RailUtils.scan(octantLoads, octantLoads, (a:Long,b:Long)=>a+b, 0L);
        
        particleData.sortAtoms(octantIds, 0, particleData.numAtoms()-1, (a:UInt,b:UInt)=> (a-b) as Int);

        val octantAtoms = new ArrayList[Pair[UInt,ArrayList[Long]]](filled);
        if (octantIds.size > 0) {
            val firstMortonId = octantIds(0);
            var currentOctant:Pair[UInt,ArrayList[Long]] = new Pair[UInt,ArrayList[Long]](firstMortonId, new ArrayList[Long](density));
            octantAtoms.add(currentOctant);

            for (i in 0..(particleData.numAtoms()-1)) {
                val leafMortonId = octantIds(i);
                if (currentOctant.first == leafMortonId) {
                    //Console.OUT.println(here + " found octant: " + leafMortonId);
                } else {
                    //Console.OUT.println(here + " creating octant: " + leafMortonId);
                    currentOctant = new Pair[UInt,ArrayList[Long]](leafMortonId, new ArrayList[Long](density));
                    octantAtoms.add(currentOctant);
                }
                //Console.OUT.println(here + " adding " + i + " to octant " + leafMortonId);
                currentOctant.second.add(i);
            }
        }

        timer.stop(FmmLocalData.TIMER_INDEX_SORT);

        return octantAtoms;
    }

    private def determineLoadBalanceLocal(octantAtoms:ArrayList[Pair[UInt,ArrayList[Long]]]) {
        val local = FastMultipoleMethod.localData;
        local.timer.start(FmmLocalData.TIMER_INDEX_BALANCE);
        val octantLoads = local.octantLoads;
        octantLoads.clear();
        for(octant in octantAtoms) {
            val leafMortonId = octant.first as Int;
            octantLoads(leafMortonId) = octant.second.size();
        }
        Team.WORLD.allreduce[Long](octantLoads, 0L, octantLoads, 0L, octantLoads.size, Team.ADD);
        val reduce = (op:(Long,Long)=>Long, init:Long, v:Rail[Long]) =>
        {
            var result:Long = init;
            for (i in 0..(v.size-1)) {
                result = op(result, v(i));
            }
            return result;
        };

        val numAtoms = reduce((a:Long,b:Long)=>a+b, 0L, octantLoads);
        val uListCost = local.cost(FmmLocalData.ESTIMATE_P2P);
        val q = numAtoms / octantLoads.size; // average particles per lowest-level octant
        val vListCost = local.cost(FmmLocalData.ESTIMATE_M2L);

        for(i in 0..(octantLoads.size-1)) {
            if (octantLoads(i) > 0) {
                val leafMortonId = i as UInt;
                val mortonId = leafMortonId | (dMax as UInt << 24);
                val id = OctantId.getFromMortonId(mortonId);
                //Console.OUT.println(id + " particles " + octantLoads(i) + " uList " + (octantLoads(i) * LeafOctant.estimateUListSize(id, dMax) * q * uListCost) + " vList size " + Octant.estimateVListSize(id, ws) + " cost " + Octant.estimateVListSize(id, ws) * vListCost);
                // octantLoad = uListInteractions * costP2P + vListInteractions * costM2L 
                octantLoads(i) = 
                        octantLoads(i) * LeafOctant.estimateUListSize(id, dMax) * q * uListCost 
                      + Octant.estimateVListSize(id, ws) * vListCost;
            }
        }

        val total = reduce((a:Long,b:Long)=>a+b, 0L, octantLoads);
        val numPlaces = Place.numPlaces();
        val placeShare = total / numPlaces;
        val leftOver = total % numPlaces;

        val firstLeafOctant = local.firstLeafOctant;
        firstLeafOctant(0) = 0UN;
        var i:Int=0n;
        for (p in 1..numPlaces) {
            val share = p <= leftOver ? placeShare+1 : placeShare;
            var load:Long = 0;
            while(load < share && i < octantLoads.size) {
                load += octantLoads(i++);
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

        local.timer.stop(FmmLocalData.TIMER_INDEX_BALANCE);
    }

    private def redistributeOctantsLocal(octantAtoms:ArrayList[Pair[UInt,ArrayList[Long]]]) {
        val local = FastMultipoleMethod.localData;
        val localOctants = local.octants;
        val particleData = local.particleData;
        local.timer.start(FmmLocalData.TIMER_INDEX_REDIST);
        val firstLeafOctant = local.firstLeafOctant;
        val levelBytes = (dMax as UInt << 24);
        finish for (var destPlace:Place=Place.places().next(here); destPlace!=here && octantAtoms.size() > 0; destPlace=Place.places().next(destPlace)) {
            val p = destPlace.id;
            val placeStart = firstLeafOctant(p);
            val placeEnd = firstLeafOctant(p+1);
            //Console.OUT.println(here + " sending octants between " + placeStart + " and " + placeEnd + " to " + destPlace);

            // redistribute octant atoms
            var end:Long = octantAtoms.size()-1;
            while (end >= 0 && octantAtoms(end).first >= placeEnd) end--;
            var start:Long = end;
            while (start >= 0 && octantAtoms(start).first >= placeStart) start--;
            ++start;

            //Console.OUT.println(here + " sending idx " + start+" ... "+end +" to " + destPlace);

            val octantsToSend = octantAtoms.moveSectionToRail(start, end);
            if (octantsToSend.size > 0) {
                //Console.OUT.println(here + " sending " + octantsToSend(0).first+" ... "+octantsToSend(octantsToSend.size-1).first+" to " + destPlace);
                val atomsStart = octantsToSend(0).second.getFirst();
                val atomsEnd = octantsToSend(octantsToSend.size-1).second.getLast();
                //Console.OUT.println(here + " sending atoms " + atomsStart + " ... " + atomsEnd + " to " + destPlace);
                val atomIds = particleData.globalIndex.moveSectionToRail(atomsStart, atomsEnd);
                val atomTypes = particleData.atomTypeIndex.moveSectionToRail(atomsStart, atomsEnd);
                val centers = particleData.x.moveSectionToRail(atomsStart, atomsEnd);
                val velocities = particleData.dx.moveSectionToRail(atomsStart, atomsEnd);
                val throwAwayForces = particleData.fx.moveSectionToRail(atomsStart, atomsEnd);

                at(destPlace) async {
                    addReceivedAtoms(atomIds, atomTypes, centers, velocities);
                }
                for (octant in octantsToSend) {
                    //Console.OUT.println(here + " removing octant " + octant.first);
                    localOctants.remove(octant.first);
                }
            }
        }
        Team.WORLD.barrier();

        // add all received particles to octants at this place
        val invSideLength = lowestLevelDim / size;
        val offset = lowestLevelDim / 2.0;
        for (i in 0..(particleData.numAtoms()-1)) {
            val leafMortonId = OctantId.getLeafMortonId(particleData.x(i), invSideLength, offset);
            val mortonId = leafMortonId | (dMax as UInt << 24);
            var leafOctant:LeafOctant = localOctants.getOrElse(mortonId, null) as LeafOctant;
            if (leafOctant == null) {
                //Console.OUT.println("octant " + OctantId.getFromMortonId(mortonId) + " stays at " + here);
                leafOctant = new LeafOctant(OctantId.getFromMortonId(mortonId), numTerms, ws, dMax);
                localOctants.put(mortonId, leafOctant);
            } else {
                //Console.OUT.println("octant " + OctantId.getFromMortonId(mortonId) + " previously created at " + here);
            }
            //Console.OUT.println(here + " adding " + particleData.globalIndex(i) + " to octant " + OctantId.getFromMortonId(mortonId));
            leafOctant.atoms.add(i);
        }
        local.timer.stop(FmmLocalData.TIMER_INDEX_REDIST);
    }

    /**
     * Add atoms received from another place to the atoms at this place.
     */
    private def addReceivedAtoms(atomIds:Rail[Long],
                                 atomTypes:Rail[Int]{self.size==atomIds.size},
                                 centers:Rail[Point3d]{self.size==atomIds.size},
                                 velocities:Rail[Vector3d]{self.size==atomIds.size}) {
        val particleData = FastMultipoleMethod.localData.particleData;
        atomic {
            for (j in 0..(atomIds.size-1)) {
                //Console.OUT.println(here + " receiving " + atomIds(j));
                particleData.addAtom(atomIds(j), atomTypes(j), centers(j), velocities(j));
            }
        }
    }

    private def createGhostChildren(parentOctant:ParentOctant) {
        val local = FastMultipoleMethod.localData;
        val octantLoads = local.octantLoads;
        val firstLeafOctant = local.firstLeafOctant;
        val id = parentOctant.id;
        var i:Long = 0;
        for (x2 in (2*id.x)..(2*id.x+1)) {
            for (y2 in (2*id.y)..(2*id.y+1)) {
                for (z2 in (2*id.z)..(2*id.z+1)) {
                    val childOctantId = OctantId(x2 as UByte, y2 as UByte, z2 as UByte, id.level+1UY);
                    val anchor = childOctantId.getAnchor(dMax);
                    var nonEmpty:Boolean = false;
                    val size = Math.pow(8.0, (dMax / childOctantId.level - 1)) as Int;
                    val startDescendant = anchor.getLeafMortonId() as Int;
                    val endDescendant = startDescendant + size - 1;
                    for (j in Math.max(firstLeafOctant(here.id+1) as Int,startDescendant)..endDescendant) {
                        if (octantLoads(j) > 0L) {
                            nonEmpty = true;
                            break;
                        }
                    }
                    if (nonEmpty) {
                        // non-empty child octant is not held at this place
                        val placeId = FastMultipoleMethod.localData.getPlaceId(anchor);
                        if (placeId != here.id) {
                            parentOctant.children(i) = new GhostOctant(childOctantId, placeId);
                       }
                    }
                    i++;
                }
            }
        }
    }

    private def createParentOctantsLocal() {
        val local = FastMultipoleMethod.localData;
        val timer = local.timer;
        timer.start(FmmLocalData.TIMER_INDEX_PARENTS);

        val octants = local.octants;
        var octantList:ArrayList[Octant] = new ArrayList[Octant](octants.size());
        for (octantEntry in octants.entries()) {
            octantList.add(octantEntry.getValue());
        }
        octantList.sort();

        // create leaf octant list
        val particleData = local.particleData;
        val leafOctants = local.leafOctants;
        leafOctants.clear();
        for (octant in octantList) {
            val leafOctant = octant as LeafOctant;
            //Console.OUT.println("at " + here + " LeafOctant " + leafOctant.id + " with " + leafOctant.atoms.size() + " atoms");
            leafOctants.add(leafOctant);
            leafOctant.makeSources(particleData);
            leafOctant.createUList(ws);
            leafOctant.createVList(ws, dMax);
        }

        // create parent octants in higher levels - assumes sorted octantList
        local.topLevelOctants.clear();
        var level:UByte=dMax;
        while (level > OctantId.TOP_LEVEL) {
            val parentOctantList = new ArrayList[Octant](octantList.size() / 8);
            var currentParent:ParentOctant = null;
            for (octant in octantList) {
                val parentId = octant.id.getParentId();
                val parentAnchorId = parentId.getAnchor(dMax);
                val placeId = local.getPlaceId(parentAnchorId);
                if (placeId == here.id) {
                    //Console.OUT.println("at " + here + " parent for " + octant.id + " = " + parentId);
                    if (currentParent == null || currentParent.id != parentId) {
                        currentParent = new ParentOctant(parentId, numTerms, ws, dMax);
                        currentParent.createVList(ws, dMax);
                        createGhostChildren(currentParent);
                        parentOctantList.add(currentParent);
                        octants.put(parentId.getMortonId(), currentParent);
                    }
                    octant.parent = currentParent;
                    val childIndex = parentId.getChildIndex(dMax, octant.id);
                    //Console.OUT.println("at " + here + " childIndex for " + octant.id + " = " + childIndex);
                    currentParent.children(childIndex) = octant;
                } else {
                    //Console.OUT.println("at " + here + " parent for " + octant.id + " = " + parentId + " anchor " + parentAnchorId + " held at " + placeId);
                    local.topLevelOctants.add(octant);
                }
            }

            octantList = parentOctantList;

            level--;
        }

        // add remaining octants to top level
        for (octant in octantList) {
            local.topLevelOctants.add(octant);
            //Console.OUT.println("at " + here + " top level " + octant.id + " " + octant.id.getMortonId());
        }

        timer.stop(FmmLocalData.TIMER_INDEX_PARENTS);
    }

    /**
     * Creates the locally essential tree at the current place.  This is later
     * used to overlap remote retrieval of particle data with other computation.
     */
    private def createLET() {
        val local = FastMultipoleMethod.localData;
        val timer = local.timer;
        timer.start(FmmLocalData.TIMER_INDEX_LET);
        val combinedUList = local.getCombinedUList(ws);
        local.locallyEssentialTree = new LET(combinedUList);
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
        val myCombinedUList = myLET.combinedUList;

        val uListPlaces = new HashMap[Long,ArrayList[UInt]](26); // a place may have up to 26 immediate neighbours
        
        // separate the uList into partial lists stored at each nearby place
        for (octantId in myCombinedUList) {
            val placeId = local.getPlaceId(octantId);
            if (placeId >= 0 && placeId < Place.numPlaces()) {
                if (placeId == here.id) {
                    // leaf octant is local
                    val octant = local.getOctant(octantId) as LeafOctant;
                    if (octant != null) {
                        myLET.setAtomDataForOctant(octantId, octant.getSources());
                    }
                } else {
                    // leaf octant exists at Place(placeId)
                    var uListForPlace : ArrayList[UInt] = uListPlaces.getOrElse(placeId, null);
                    if (uListForPlace == null) {
                        uListForPlace = new ArrayList[UInt]();
                        uListPlaces.put(placeId, uListForPlace);
                    }
                    uListForPlace.add(octantId);
                }
            }
        }

        // retrieve the partial list for each remote place and store into my LET
        finish for (placeEntry in uListPlaces.entries()) async {
            val placeId = placeEntry.getKey();
            val uListForPlace = placeEntry.getValue();
            val uListArray = uListForPlace.toRail();
            val atomsForPlace = at(Place(placeId)) { FastMultipoleMethod.getAtomsForOctantList(uListArray)};
            for (i in 0..(uListArray.size-1)) {
                myLET.setAtomDataForOctant(uListArray(i), atomsForPlace(i));
            }
        }
        local.timer.stop(FmmLocalData.TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of octant indexes stored at a single
     * place, returns an Rail, each element of which is in turn
     * a Rail[Double] containing the atoms for each octant.
     */
    private static def getAtomsForOctantList(octantList:Rail[UInt]) {
        val local = FastMultipoleMethod.localData;
        val atomList = new Rail[Rail[Double]](octantList.size);
        for (i in 0..(octantList.size-1)) {
            val octant = local.getOctant(octantList(i)) as LeafOctant;
            if (octant != null) {
                atomList(i) = octant.getSources();
            }
        }
        return atomList;
    }
}

