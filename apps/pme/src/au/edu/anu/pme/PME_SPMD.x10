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
package au.edu.anu.pme;

import x10.compiler.Inline;
import x10.regionarray.Array;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;
import x10.regionarray.PeriodicDist;
import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.Pair;
import x10.util.Team;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.PointCharge;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.fft.DistributedReal3dFft_SPMD;
import au.edu.anu.util.Timer;
//import org.netlib.fdlibm.Erf;

/**
 * This class implements a Smooth Particle Mesh Ewald method to calculate
 * the potential of a system of charged particles.
 * TODO dipole correction for forces as described in Lambert, Darden & Board and Deem et al.
 * @see Essmann et al. "A Smooth Particle Mesh Ewald method", J. Comp. Phys. 101,
 * pp.8577-8593 (1995) DOI: 10.1063/1.470117
 */
public class PME_SPMD {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int           = 0n;
    public static val TIMER_INDEX_GRIDCHARGES : Int     = 1n;
    public static val TIMER_INDEX_INVFFT : Int          = 2n;
    public static val TIMER_INDEX_THETARECCONVQ : Int   = 3n;
    public static val TIMER_INDEX_RECIPROCAL : Int      = 4n;
    public static val TIMER_INDEX_SELF : Int            = 5n;
    public static val TIMER_INDEX_DIRECT : Int          = 6n;
    public static val TIMER_INDEX_PREFETCH : Int        = 7n;
    public static val TIMER_INDEX_SETUP : Int           = 8n;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    static val timer = new Timer(9);

    /** The number of grid lines in each dimension of the simulation unit cell. */
    private val gridSize : Rail[Long];

    /** Double representations of the various grid dimensions */
    private val K1 : Double;
    private val K2 : Double; 
    private val K3 : Double;

    /** The edges of the unit cell. */
    private val edges : Rail[Vector3d];
    private val edgeLengths : Rail[Double];

    /** The conjugate reciprocal vectors for each dimension. */
    private val edgeReciprocals : Rail[Vector3d];

    /** 
     * Scaling vector for scaled fractional coordinates.
     * Assumes unit rectangular cell. 
     */
    private val scalingVector : Vector3d;

    private val gridDist : Dist(3);
    
    /** The order of B-spline interpolation */
    private val splineOrder : Int;

    /** The Ewald coefficient beta */
    private val beta : Double;

    /** The direct sum cutoff distance in Angstroms */
    private val cutoff : Double;

    /** 
     * Translation vectors for neighbouring unit cells 
     * (the 26 cells surrounding the origin cell).
     * These are replicated across all places using a PlaceLocalHandle.
     * Dimensions of the enclosed array are:
     * 0: x translation (difference between x-coordinate of sub-cells
     * 1: y translation
     * 2: z translation
     */
    private val imageTranslations : PlaceLocalHandle[Array[Vector3d](3){rect}];

    /** The atoms in the simulation, divided up into a distributed array of Arrays, one for each place. */
    private val atoms : DistArray[Rail[MMAtom]](1);

    private val B : DistArray[Double]{self.dist==gridDist};
    private val C : DistArray[Double]{self.dist==gridDist};
    private val BdotC : DistArray[Double]{self.dist==gridDist};

    /** The gridded charge array Q as defined in Eq. 4.6 */
    private val Q : DistArray[Double]{self.dist==gridDist};

    /** The inverse DFT of the Q array.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val Qinv : DistArray[Complex]{self.dist==gridDist};

    /** thetaRecConvQ as used in Eq. 4.7.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val thetaRecConvQ : DistArray[Complex]{self.dist==gridDist};
    private val thetaRecConvQReal : DistArray[Double]{self.dist==gridDist};

    /** Scratch array for use during 3D FFT.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val temp : DistArray[Complex]{self.dist==gridDist};
    private val temp2 : DistArray[Double]{self.dist==gridDist};

    /** 
     * An array of box divisions within the unit cell, with a side length
     * equal to half the direct sum cutoff distance.  (N.B. if the unit cell
     * side length is not an exact multiple of the subcell side length, the
     * last box in each dimension will be smaller than the cutoff distance, 
     * resulting in anisotropy in the direct potential.)
     * Direct sums are only calculated between particles in the same box and
     * the 26 neighbouring boxes.
     * Dimensions of the array region are (x,y,z)
     * TODO assumes cubic unit cell
     */
    private val subCells : DistArray[Rail[PointCharge]](3);
    /** The number of sub-cells per side of the unit cell. */
    private val numSubCells : Long;

    /** 
     * A cache of atoms from subcells stored at other places.  
     * This is used to prefetch atom data for direct energy calculation.
     */
    private val atomsCache : DistArray[Array[Rail[PointCharge]]{rank==3,rect}](1);

    /**
     * Creates a new particle mesh Ewald method.
     * @param edges the edge vectors of the unit cell
     * @param gridSize the number of grid lines in each dimension of the unit cell
     * @param atoms the atoms in the unit cell
     * @param splineOrder the order n of B-spline interpolationb
     * @param beta the Ewald coefficient beta
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(edges : Rail[Vector3d],
            gridSize : Rail[Long],
            atoms: DistArray[Rail[MMAtom]](1),
            splineOrder : Int,
            beta : Double,
            cutoff : Double) {
        this.gridSize = gridSize;
        val K1 = gridSize(0) as Double;
        val K2 = gridSize(1) as Double;
        val K3 = gridSize(2) as Double;
        this.edges = edges;
        this.edgeLengths = new Rail[Double](3, (i:Long) => edges(i).length());
        this.edgeReciprocals = new Rail[Vector3d](3, (i:Long) => edges(i).inverse());
        this.scalingVector = Vector3d(edges(0).inverse().i * K1, edges(1).inverse().j * K2, edges(2).inverse().k * K3);
        this.K1 = K1;
        this.K2 = K2;
        this.K3 = K3;

        this.atoms = atoms;
        val gridRegion = Region.make(0..(gridSize(0)-1), 0..(gridSize(1)-1), 0..(gridSize(2)-1));
        val gridDist = Dist.makeBlockBlock(gridRegion, 0, 1);
        this.gridDist = gridDist;
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        this.imageTranslations = PlaceLocalHandle.make[Array[Vector3d](3){rect}](
            PlaceGroup.WORLD, 
            () => new Array[Vector3d](Region.make(-1..1, -1..1, -1..1), 
                ([i,j,k] : Point(3)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k))) 
        );

        if (edgeLengths(0) % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + edgeLengths(0) + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(edgeLengths(0) / (cutoff/2.0)) as Long;
        val subCellRegion = Region.make(0..(numSubCells-1), 0..(numSubCells-1), 0..(numSubCells-1));
        val subCells = DistArray.make[Rail[PointCharge]](new PeriodicDist(Dist.makeBlockBlock(subCellRegion, 0, 1)));
        //Console.OUT.println("subCells dist = " + subCells.dist);
        this.subCells = subCells;
        this.numSubCells = numSubCells;

        val atomsCache = DistArray.make[Array[Rail[PointCharge]]{rank==3,rect}](Dist.makeUnique());
        finish ateach(p in atomsCache) {
            val mySubCellRegion = subCells.dist(here);
            if (! mySubCellRegion.isEmpty()) {
                val directRequiredRegion = 
                    Region.make((mySubCellRegion.min(0) - 2)..(mySubCellRegion.max(0) + 2),
                                (mySubCellRegion.min(1) - 2)..(mySubCellRegion.max(1) + 2),
                                (mySubCellRegion.min(2) - 2)..(mySubCellRegion.max(2) + 2));
                atomsCache(p) = new Array[Rail[PointCharge]](directRequiredRegion);
            }
        }
        this.atomsCache = atomsCache;

        //Console.OUT.println("gridDist = " + gridDist);

        Q = DistArray.make[Double](gridDist);
        BdotC = DistArray.make[Double](gridDist);

        // TODO following arrays should be scoped local variables XTENLANG-???
        Qinv = DistArray.make[Complex](gridDist);
        thetaRecConvQ = DistArray.make[Complex](gridDist);
        thetaRecConvQReal = DistArray.make[Double](gridDist);
        temp = DistArray.make[Complex](gridDist);
        temp2 = DistArray.make[Double](gridDist);
        B = DistArray.make[Double](gridDist);
        C = DistArray.make[Double](gridDist);
    }

    /**
     * This method sets up the B, C and BdotC arrays, which only need be done once
     * in the simulation (as it is dependent on the grid size).
     * TODO this should be done in the constructor, but can't be done without big 
     * memory leaks due to proto rules.
     */
    public def setup() {
        timer.start(TIMER_INDEX_SETUP);
        initBArray();
        initCArray();
        B.map(BdotC, C, (a:Double, b:Double) => (a*b));
        divideAtomsIntoSubCells();
        timer.stop(TIMER_INDEX_SETUP);
    }
	
    public def getEnergy() : Double {
        val totalEnergy = finish (Reducible.SumReducer[Double]()) {
            ateach(p1 in Dist.makeUnique()) {
                timer.start(TIMER_INDEX_TOTAL);
                prefetchPackedAtomsLocal();

                gridChargesLocal();
                Team.WORLD.barrier();

                timer.start(TIMER_INDEX_INVFFT);
                DistributedReal3dFft_SPMD.doFFT3d(Q, Qinv, temp);
                Team.WORLD.barrier();
                timer.stop(TIMER_INDEX_INVFFT);

                timer.start(TIMER_INDEX_THETARECCONVQ);
                // create F^-1(thetaRecConvQ)
                val localBdotC = BdotC.getLocalPortion();
                val localThetaRecConvQ = thetaRecConvQ.getLocalPortion();
                val localQinv = Qinv.getLocalPortion();
                val localRegion = localBdotC.region as Region(3){rect};
                localBdotC.map[Complex,Complex](localThetaRecConvQ, localQinv, (a:Double, b:Complex) => a * b);

                Team.WORLD.barrier();
                // and do inverse FFT
                DistributedReal3dFft_SPMD.doFFT3d(thetaRecConvQ, thetaRecConvQReal, temp, temp2);
                timer.stop(TIMER_INDEX_THETARECCONVQ);

                val reciprocalEnergy = getReciprocalEnergy(thetaRecConvQReal);
                val selfEnergy = getSelfEnergy();
                val directEnergy = getDirectEnergy();
                //Console.OUT.println("directEnergy = " + directEnergy);
                //Console.OUT.println("selfEnergy = " + selfEnergy);
                //Console.OUT.println("correctionEnergy = " + correctionEnergy);
                val correctionEnergy = 0.0;
                //Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
                
                val myEnergy = directEnergy + reciprocalEnergy + (correctionEnergy + selfEnergy);

                timer.stop(TIMER_INDEX_TOTAL);
                Team.WORLD.allreduce[Long](timer.total, 0L, timer.total, 0L, timer.total.size, Team.MAX);

                offer myEnergy;
            }
        };


        return totalEnergy;
    }

    /**
     * Divide the atoms into a grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def divideAtomsIntoSubCells() {
        val halfCutoff = (cutoff / 2.0);
        val subCellsTemp = DistArray.make[ArrayList[PointCharge]](subCells.dist, (Point) => new ArrayList[PointCharge]());
        val atoms = this.atoms; // TODO shouldn't be necessary XTENLANG-1913
        val halfNumSubCells = this.numSubCells / 2;
        finish ateach(p in atoms) {
            val localAtoms = atoms(p);
            finish for (l in 0..(localAtoms.size-1)) {
                val atom = localAtoms(l);
                val centre = atom.centre;
                val charge = atom.charge;
                // get subcell i,j,k
                val i = (centre.i / halfCutoff) as Int + halfNumSubCells;
                val j = (centre.j / halfCutoff) as Int + halfNumSubCells;
                val k = (centre.k / halfCutoff) as Int + halfNumSubCells;
                at(subCellsTemp.dist(i,j,k)) async {
                    atomic subCellsTemp(i,j,k).add(new PointCharge(centre, charge));
                }
            }
        }
        val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach([i,j,k] in subCells) {
            subCells(i,j,k) = subCellsTemp(i,j,k).toRail();
        }
    }

    /**
     * At each place, fetch all required atoms from neighbouring
     * places for direct calculation.
     */
    private def prefetchPackedAtomsLocal() {
        timer.start(TIMER_INDEX_PREFETCH);
        val myAtomsCache = atomsCache(here.id);
        if (myAtomsCache != null) {
            val haloPlaces = new HashMap[Long,ArrayList[Point(3)]](8); // a place may have up to 8 immediate neighbours in the two block-divided dimensions
            
            // separate the halo subcells into partial lists stored at each nearby place
            for (boxIndex[x,y,z] in myAtomsCache.region) {
                val placeId = subCells.dist(x,y,z).id;
                var haloForPlace : ArrayList[Point(3)] = haloPlaces.getOrElse(placeId, null);
                if (haloForPlace == null) {
                    haloForPlace = new ArrayList[Point(3)]();
                    haloPlaces.put(placeId, haloForPlace);
                }
                haloForPlace.add(boxIndex);
            }

            // retrieve the partial list for each place and store into my LET
            finish for (placeEntry in haloPlaces.entries()) async {
                val placeId = placeEntry.getKey();
                val haloForPlace = placeEntry.getValue();
                val haloListArray = haloForPlace.toRail();
                if (placeId == here.id) {
                    // atoms cache is just a set of pointers to sub cells that are here
                    for (i in 0..(haloListArray.size-1)) {
                        myAtomsCache(haloListArray(i)) = subCells(haloListArray(i));
                    }
                } else {
                    val atomsForPlace = at(Place.place(placeId)) { getAtomsForSubcellList(subCells, haloListArray)};
                    for (i in 0..(haloListArray.size-1)) {
                        myAtomsCache(haloListArray(i)) = atomsForPlace(i);
                    }
                }
            }
        }
        timer.stop(TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of subcell indices as Point(3) stored at a single
     * place, returns an Array, each element of which is in turn
     * a Array of PointCharge containing the atoms for each subcell.
	 * TODO should use instance fields instead of all those parameters - XTENLANG-1913
     */
    private static def getAtomsForSubcellList(subCells : DistArray[Rail[PointCharge]](3), boxList : Rail[Point(3)]) {
        val atoms = new Rail[Rail[PointCharge]](boxList.size, 
                                                 (i:Long) => subCells(boxList(i)));
        return atoms;
    }

    public def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        val cutoffSquared = cutoff*cutoff;
		val subCellsDist = this.subCells.dist;
        val cachedAtoms = atomsCache(here.id);
        val translations = imageTranslations();
        val localRegion = subCellsDist(here) as Region(3){rect};
        val directEnergy = finish (Reducible.SumReducer[Double]()) {
            for ([x,y,z] in localRegion) async {
                val thisCell = cachedAtoms(x,y,z) as Rail[PointCharge];
                var cellDirectEnergy : Double = 0.0;
                for (var i:Long = x-2; i<=x; i++) {
                    var n1:Int = 0n;
                    if (i < 0) {
                        n1 = -1n;
                    } // can't have (i > numSubCells+1)
                    for (var j:Long = y-2; j<=y+2; j++) {
                        var n2:Int = 0n;
                        if (j < 0) {
                            n2 = -1n;
                        } else if (j > numSubCells-1) {
                            n2 = 1n;
                        }
                        for (var k:Long = z-2; k<=z+2; k++) {
                            var n3:Int = 0n;
                            if (k < 0) {
                                n3 = -1n;
                            } else if (k > numSubCells-1) {
                                n3 = 1n;
                            }
                            // interact with "left half" of other boxes i.e. only boxes with i<=x
                            if (i < x || (i == x && j < y) || (i == x && j == y && k < z)) {
                                val translation = translations(n1,n2,n3);
                                val otherCell : Rail[PointCharge] = cachedAtoms(i,j,k);
                                for (otherAtomIndex in 0..(otherCell.size-1)) {
                                    val otherAtom = otherCell(otherAtomIndex);
                                    val imageLoc = otherAtom.centre + translation;
                                    val otherAtomCharge = otherAtom.charge;
                                    for (thisAtomIndex in 0..(thisCell.size-1)) {
                                        val thisAtom = thisCell(thisAtomIndex);
                                        val rSquared = thisAtom.centre.distanceSquared(imageLoc);
                                        if (rSquared < cutoffSquared) {
                                            val r = Math.sqrt(rSquared);
                                            val chargeProduct = thisAtom.charge * otherAtomCharge;
                                            val imageDirectComponent = chargeProduct * Math.erfc(beta * r) / r;
                                            cellDirectEnergy += imageDirectComponent;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // atoms in same cell
                for (i in 0..(thisCell.size-1)) {
                    val thisAtom = thisCell(i);
                    for (j in 0..(i-1)) {
                        val otherAtom = thisCell(j);
                        val rjri = otherAtom.centre - thisAtom.centre;
                        val rSquared = rjri.lengthSquared();
                        if (rSquared < cutoffSquared) {
                            val r = Math.sqrt(rSquared);
                            val directComponent = thisAtom.charge * otherAtom.charge * Math.erfc(beta * r) / r;
                            cellDirectEnergy += directComponent;
                        }
                    }
                }
                offer cellDirectEnergy;
            }
        };
        timer.stop(TIMER_INDEX_DIRECT);
        return directEnergy;
    }

    /**
     * Returns the self energy as defined in Eq. 2.5, which is the
     * contribution of the interaction in reciprocal space of each
     * atom with itself.  This is subtracted from the final energy.
     */
    public def getSelfEnergy() : Double {
        timer.start(TIMER_INDEX_SELF);
        val localSubCells = subCells.getLocalPortion();
        val localRegion = localSubCells.region as Region(3){rect};
        var mySelfEnergy : Double = 0.0;
        for ([i,j,k] in localRegion) {
            val thisCell = localSubCells(i,j,k) as Rail[PointCharge];
            for (thisAtom in 0..(thisCell.size-1)) {
                mySelfEnergy += thisCell(thisAtom).charge * thisCell(thisAtom).charge;
            }
        }
        val selfEnergy = mySelfEnergy * -beta / Math.sqrt(Math.PI);
        timer.stop(TIMER_INDEX_SELF);
        return selfEnergy;
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     * Spline interpolation code inspired by GROMACS http://www.gromacs.org/
     */
    public def gridChargesLocal() {
        timer.start(TIMER_INDEX_GRIDCHARGES);

        val qLocal = Q.getLocalPortion() as Array[Double](3){rect};
        val localGridRegion = qLocal.region as Region(3){rect};
        if (!localGridRegion.isEmpty()) {
            val subCellRegion = subCells.dist.region;
            qLocal.clear();
            val gridSize0 = gridSize(0);
            val gridSize1 = gridSize(1);
            val gridSize2 = gridSize(2);
            val myAtomsCache = atomsCache(here.id);

            val subCellHaloRegion = PME_SPMD.getSubcellHaloRegionForPlace(gridSize0, numSubCells, splineOrder, localGridRegion, subCellRegion);
            //Console.OUT.println("subCellHaloRegion at " + here + " = " + subCellHaloRegion);
            val iSpline = new Rail[Double](splineOrder);
            val jSpline = new Rail[Double](splineOrder);
            val kSpline = new Rail[Double](splineOrder);

            val qiMin = localGridRegion.min(0);
            val qiMax = localGridRegion.max(0);
            val qjMin = localGridRegion.min(1);
            val qjMax = localGridRegion.max(1);

            for ([x,y,z] in subCellHaloRegion) {
                val thisCell = myAtomsCache(x,y,z) as Rail[PointCharge];
                for (atomIndex in 0..(thisCell.size-1)) {
                    val atom = thisCell(atomIndex);
                    val q = atom.charge;
                    val u = atom.centre.scale(scalingVector) + Vector3d(K1, K2, K3); // Eq. 3.1 TODO general non-cubic
                    val u1i = u.i as Int;
                    PME_SPMD.fillSpline(1.0 - (u.i - u1i), iSpline, splineOrder);
                    val u2i = u.j as Int;
                    PME_SPMD.fillSpline(1.0 - (u.j - u2i), jSpline, splineOrder);
                    val u3i = u.k as Int;
                    PME_SPMD.fillSpline(1.0 - (u.k - u3i), kSpline, splineOrder);
                    for (i in 0..(splineOrder-1)) {
                        val k1 = (u1i - i + gridSize0) % gridSize0;
                        val iVal = q * iSpline(i);
                        if (k1 >= qiMin && k1 <= qiMax) {
                            for (j in 0..(splineOrder-1)) {
                                val k2 = (u2i - j + gridSize1) % gridSize1;
                                val jVal = iVal * jSpline(j);
                                if (k2 >= qjMin && k2 <= qjMax) {
                                    for (k in 0..(splineOrder-1)) {
                                        val k3 = (u3i - k + gridSize2) % gridSize2;
                                        val kVal = jVal * kSpline(k);
                                        qLocal(k1,k2,k3) += kVal;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        timer.stop(TIMER_INDEX_GRIDCHARGES);
    }

    /**
     * Given the grid region for a particular place, returns the Region of
     * subCells that will contribute to the grid at that place.
     * Some of this region will be resident at the given place; other parts
     * may be held as ghost cells.
     */
    private static def getSubcellHaloRegionForPlace(gridSize:Long, numSubCells:Long, splineOrder:Int, gridRegion:Region(3){rect}, subCellRegion:Region(3){rect}) {
        val gridPointsPerSubCell = (gridSize as Double) / numSubCells;

        val subCellMinX = Math.floor(gridRegion.min(0) / gridPointsPerSubCell) as Int;

        val fullX = gridRegion.min(0) == 0L && gridRegion.max(0) == (gridSize-1);
        val subCellMaxX:Long;
        if (fullX) {
            subCellMaxX = subCellRegion.max(0);
        } else {
            val maxX = gridRegion.max(0) + (splineOrder-1);
            subCellMaxX = Math.ceil(maxX / gridPointsPerSubCell) as Long;
        }

        val subCellMinY = Math.floor(gridRegion.min(1) / gridPointsPerSubCell) as Long;
        val fullY = gridRegion.min(1) == 0L && gridRegion.max(1) == (gridSize-1);
        val subCellMaxY:Long;
        if (fullY) {
            subCellMaxY = subCellRegion.max(1);
        } else {
            val maxY = gridRegion.max(1) + (splineOrder-1);
            subCellMaxY = Math.ceil(maxY / gridPointsPerSubCell) as Long;
        }
        // each place holds a full pencil of subcells through Z dimension
        val subCellHaloRegion = 
            Region.make((subCellMinX..subCellMaxX),
                        (subCellMinY..subCellMaxY),
                        (subCellRegion.min(2)..subCellRegion.max(2)));
        return subCellHaloRegion;
    }

    /**
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(thetaRecConvQ:DistArray[Double]{self.dist==gridDist}) {
        timer.start(TIMER_INDEX_RECIPROCAL);

        val scale = 1.0 / (K1 * K2 * K3);
        var myReciprocalEnergy : Double = 0.0;
        val localQ = Q.getLocalPortion();
        val localThetaRecConvQ = thetaRecConvQ.getLocalPortion();
        val localRegion = localQ.region as Region(3){rect};
        for ([i,j,k] in localRegion) {
            myReciprocalEnergy += localQ(i,j,k) * localThetaRecConvQ(i,j,k);
        }

        val reciprocalEnergy = scale * myReciprocalEnergy / 2.0;
        timer.stop(TIMER_INDEX_RECIPROCAL);

        return reciprocalEnergy;
    }

    private static @Inline def fillSpline(offset:Double, spline:Rail[Double], splineOrder:Int) {
        spline(splineOrder-1) = 0.0;
        spline(1) = offset;
        spline(0) = 1.0 - offset;
    
        for(k in 3..(splineOrder-1)) {
            val div = 1.0 / (k-1.0);    
            spline(k-1) = div * offset * spline(k-2);
            for(l in 1..(k-2)) {
                spline(k-l-1) = div * ((offset+l) * spline(k-l-2) + (k-l-offset) * spline(k-l-1));
            }
            spline(0) = div * (1.0-offset) * spline(0);
        }
    
        val div = 1.0 / (splineOrder-1);
        spline(splineOrder-1) = div * offset * spline(splineOrder-2);
        for(l in 1..(splineOrder-2)) {
            spline(splineOrder-l-1) = div * 
                ((offset+l) * spline(splineOrder-l-2) + (splineOrder-l-offset) * spline(splineOrder-l-1));
        }
        spline(0) = div * (1.0-offset) * spline(0);
    }

    /**
     * Initialises the array B as defined by Eq 4.8 and 4.4
     * TODO should be able to construct array as local and return, but can't due to GC limitations
     */
    public def initBArray() {
        val B = this.B; // TODO shouldn't be necessary XTENLANG-1913
        val splineOrder = this.splineOrder; // TODO shouldn't be necessary XTENLANG-1913
		val K1 = this.K1; // TODO shouldn't be necessary XTENLANG-1913;
		val K2 = this.K2; // TODO shouldn't be necessary XTENLANG-1913;
		val K3 = this.K3; // TODO shouldn't be necessary XTENLANG-1913;
        finish ateach(place1 in Dist.makeUnique()) {
            val splines = new Array[Double](splineOrder);
            for (k in 1..(splineOrder-1)) {
                splines(k) = bSpline(splineOrder, k);
            }
            val regionHere = B.dist(here) as Region(3){rect};
            for ([m1,m2,m3] in regionHere) {
                val m1D = m1 as Double;
                var sumK1 : Complex = Complex.ZERO;
                for (k in 0..(splineOrder-2)) {
                    sumK1 = sumK1 + splines(k+1) * Math.exp(2.0 * Math.PI * m1D * k / K1 * Complex.I);
                }
                val b1 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m1D / K1 * Complex.I) / sumK1).abs();

                val m2D = m2 as Double;
                var sumK2 : Complex = Complex.ZERO;
                for (k in 0..(splineOrder-2)) {
                    sumK2 = sumK2 + splines(k+1) * Math.exp(2.0 * Math.PI * m2D * k / K2 * Complex.I);
                }
                val b2 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m2D / K2 * Complex.I) / sumK2).abs();
                    
                val m3D = m3 as Double;
                var sumK3 : Complex = Complex.ZERO;
                for (k in 0..(splineOrder-2)) {
                    sumK3 = sumK3 + splines(k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
                }
                val b3 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m3D / K3 * Complex.I) / sumK3).abs();

                B(m1,m2,m3) = b1 * b1 * b2 * b2 * b3 * b3;
                //Console.OUT.println("B(" + m1 + "," + m2 + "," + m3 + ") = " + B(m1,m2,m3));
            }
        }
    }

    /**
     * Initialises the array C as defined by Eq 3.9
     * TODO should be able to construct array as local and return, but can't due to GC limitations
     */
    public def initCArray() {
        val V = getVolume();
        //Console.OUT.println("V = " + V);
        val C = this.C; // TODO shouldn't be necessary XTENLANG-1913
		val K1 = this.K1; // TODO shouldn't be necessary XTENLANG-1913;
		val K2 = this.K2; // TODO shouldn't be necessary XTENLANG-1913;
		val K3 = this.K3; // TODO shouldn't be necessary XTENLANG-1913;
        val edgeReciprocals = this.edgeReciprocals; // TODO shouldn't be necessary XTENLANG-1913
        val beta = this.beta; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach(place1 in Dist.makeUnique()) {
            val regionHere = C.dist(here) as Region(3){rect};
            for ([m1,m2,m3] in regionHere) {
                val m1prime = m1 <= K1/2 ? m1 as Double : m1 - K1;
                val m2prime = m2 <= K2/2 ? m2 as Double : m2 - K2;
                val m3prime = m3 <= K3/2 ? m3 as Double : m3 - K3;
                val mVec = edgeReciprocals(0).mul(m1prime).add(edgeReciprocals(1).mul(m2prime)).add(edgeReciprocals(2).mul(m3prime));
                val mSquared = mVec.dot(mVec);
                C(m1,m2,m3) = Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)) / (mSquared * Math.PI * V);
            }
        }
        at(C.dist(0,0,0)) {
            C(0,0,0) = 0.0;
        }
    }

    /* 
     * Gets the nth order B-spline M_n(u) as per Eq. 4.1
     */
    public final static def bSpline(n : Int, u : Double) : Double {
        if (n == 4n) {
            return bSpline4(u);
        } else if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 2n) {
            return 1.0 - Math.abs(u - 1.0);
        } else {
            return u / (n-1n) * bSpline(n-1n, u) + (n - u) / (n - 1n) * bSpline(n-1n, u-1.0);
        }
    }

    /* 
     * Gets the 4th order B-spline M_4(u) as per Eq. 4.1
     */
    private final static def bSpline4(u : Double) : Double {
        if (u <= 0.0 || u >= 4) {
            return 0.0;
        } else {
            return u / 3 * bSpline3(u) + (4 - u) / 3 * bSpline3(u-1.0);
        }
    }

    /* 
     * Gets the 3rd order B-spline M_3(u) as per Eq. 4.1
     */
    private final static def bSpline3(u : Double) : Double {
        if (u <= 0.0 || u >= 3) {
            return 0.0;
        } else {
            return u / 2 * bSpline2(u) + (3 - u) / 2 * bSpline2(u-1.0);
        }
    }

    /* 
     * Gets the 2nd order B-spline M_2(u) as per Eq. 4.1
     */
    private final static def bSpline2(u : Double) : Double {
        if (u <= 0.0 || u >= 2) {
            return 0.0;
        } else {
            return 1.0 - Math.abs(u - 1.0);
        }
    }

    /**
     * Gets the volume V of the unit cell.
     */
    private def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }
}

