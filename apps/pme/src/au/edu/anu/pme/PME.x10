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

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.PointCharge;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.fft.Distributed3dFft;
import au.edu.anu.util.Timer;
//import org.netlib.fdlibm.Erf;

import x10.compiler.Inline;
import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.Pair;

/**
 * This class implements a Smooth Particle Mesh Ewald method to calculate
 * the potential of a system of charged particles.
 * TODO dipole correction for forces as described in Lambert, Darden & Board and Deem et al.
 * @see Essmann et al. "A Smooth Particle Mesh Ewald method", J. Comp. Phys. 101,
 * pp.8577-8593 (1995) DOI: 10.1063/1.470117
 */
public class PME {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int = 0;
    public static val TIMER_INDEX_PREFETCH : Int = 1;
    public static val TIMER_INDEX_DIRECT : Int = 2;
    public static val TIMER_INDEX_SELF : Int = 3;
    public static val TIMER_INDEX_GRIDCHARGES : Int = 4;
    public static val TIMER_INDEX_INVFFT : Int = 5;
    public static val TIMER_INDEX_THETARECCONVQ : Int = 6;
    public static val TIMER_INDEX_RECIPROCAL : Int = 7;
    public static val TIMER_INDEX_SETUP : Int = 8;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(9);

    /** The number of grid lines in each dimension of the simulation unit cell. */
    private val gridSize : Rail[Int];

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
     * These are replicated across all places in a DistArray with a unique dist.
     * Dimensions of the enclosed array are:
     * 0: x translation (difference between x-coordinate of sub-cells
     * 1: y translation
     * 2: z translation
     */
    private val imageTranslations : DistArray[Array[Vector3d](3){rect}](1);

    /** The atoms in the simulation, divided up into a distributed array of Arrays, one for each place. */
    private val atoms : DistArray[Rail[MMAtom]](1);

    private val B : DistArray[Double]{self.dist==gridDist};
    private val C : DistArray[Double]{self.dist==gridDist};
    private val BdotC : DistArray[Double]{self.dist==gridDist};

    /** The gridded charge array Q as defined in Eq. 4.6 */
    private val Q : DistArray[Complex]{self.dist==gridDist};

    /** The inverse DFT of the Q array.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val Qinv : DistArray[Complex]{self.dist==gridDist};

    /** thetaRecConvQ as used in Eq. 4.7.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val thetaRecConvQ : DistArray[Complex]{self.dist==gridDist};

    /** Scratch array for use during 3D FFT.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val temp : DistArray[Complex]{self.dist==gridDist};

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
    private val numSubCells : Int;

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
            gridSize : Rail[Int],
            atoms: DistArray[Rail[MMAtom]](1),
            splineOrder : Int,
            beta : Double,
            cutoff : Double) {
        this.gridSize = gridSize;
        val K1 = gridSize(0) as Double;
        val K2 = gridSize(1) as Double;
        val K3 = gridSize(2) as Double;
        this.edges = edges;
        this.edgeLengths = new Array[Double](3, (i : Int) => edges(i).length());
        this.edgeReciprocals = new Array[Vector3d](3, (i : Int) => edges(i).inverse());
        this.scalingVector = Vector3d(edges(0).inverse().i * K1, edges(1).inverse().j * K2, edges(2).inverse().k * K3);
        this.K1 = K1;
        this.K2 = K2;
        this.K3 = K3;

        this.atoms = atoms;
        val gridRegion = (0..(gridSize(0)-1)) * (0..(gridSize(1)-1)) * (0..(gridSize(2)-1));
        gridDist = Dist.makeBlockBlock(gridRegion, 0, 1);
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        this.imageTranslations = DistArray.make[Array[Vector3d](3){rect}](
            Dist.makeUnique(), 
            new Array[Vector3d](-1..1 * -1..1 * -1..1, 
                ([i,j,k] : Point(3)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k))) 
        );

        if (edgeLengths(0) % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + edgeLengths(0) + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(edgeLengths(0) / (cutoff/2.0)) as Int;
        val subCellRegion = 0..(numSubCells-1) * 0..(numSubCells-1) * 0..(numSubCells-1);
        val subCells = DistArray.make[Rail[PointCharge]](new PeriodicDist(Dist.makeBlockBlock(subCellRegion, 0, 1)));
        //Console.OUT.println("subCells dist = " + subCells.dist);
        this.subCells = subCells;
        this.numSubCells = numSubCells;

        val atomsCache = DistArray.make[Array[Rail[PointCharge]]{rank==3,rect}](Dist.makeUnique());
        finish ateach(p in atomsCache) {
            val mySubCellRegion = (subCells.dist | here).region;
            if (! mySubCellRegion.isEmpty()) {
                val directRequiredRegion = ((mySubCellRegion.min(0) - 2)..(mySubCellRegion.max(0) + 2))
                                         * ((mySubCellRegion.min(1) - 2)..(mySubCellRegion.max(1) + 2))
                                         * ((mySubCellRegion.min(2) - 2)..(mySubCellRegion.max(2) + 2));
                atomsCache(p) = new Array[Rail[PointCharge]](directRequiredRegion);
            }
        }
        this.atomsCache = atomsCache;

        //Console.OUT.println("gridDist = " + gridDist);

        Q = DistArray.make[Complex](gridDist);
        BdotC = DistArray.make[Double](gridDist);

        // TODO following arrays should be scoped local variables XTENLANG-???
        Qinv = DistArray.make[Complex](gridDist);
        thetaRecConvQ = DistArray.make[Complex](gridDist);
        temp = DistArray.make[Complex](gridDist);
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
        timer.start(TIMER_INDEX_TOTAL);

        finish {
            async { prefetchRemoteAtoms(); }

            gridCharges();

            timer.start(TIMER_INDEX_INVFFT);
            new Distributed3dFft(gridSize(0), Q, Qinv, temp).doFFT3d(false);
            timer.stop(TIMER_INDEX_INVFFT);

            timer.start(TIMER_INDEX_THETARECCONVQ);
            // create F^-1(thetaRecConvQ)
            // TODO DistArray.map uses slow apply,set operators
            //BdotC.map(thetaRecConvQ, Qinv, (a:Double, b:Complex) => (a*b));
            val BdotC = this.BdotC; // TODO shouldn't be necessary XTENLANG-1913
            val thetaRecConvQ = this.thetaRecConvQ; // TODO shouldn't be necessary XTENLANG-1913
            val Qinv = this.Qinv; // TODO shouldn't be necessary XTENLANG-1913
            finish ateach(placeId in Dist.makeUnique(gridDist.places())) {
                val localBdotC = BdotC.getLocalPortion();
                val localThetaRecConvQ = thetaRecConvQ.getLocalPortion();
                val localQinv = Qinv.getLocalPortion();
                val localRegion = localBdotC.region as Region(3){rect};
                for ([x,y,z] in localRegion) {
                    localThetaRecConvQ(x,y,z) = localBdotC(x,y,z) * localQinv(x,y,z);
                }
            }

            // and do inverse FFT
            new Distributed3dFft(gridSize(0), thetaRecConvQ, thetaRecConvQ, temp).doFFT3d(true);
            timer.stop(TIMER_INDEX_THETARECCONVQ);
        }

        val reciprocalEnergy = getReciprocalEnergy(thetaRecConvQ);
        val selfEnergy = getSelfEnergy();
        val directEnergy = getDirectEnergy();

        //Console.OUT.println("directEnergy = " + directEnergy);
        //Console.OUT.println("selfEnergy = " + selfEnergy);
        //Console.OUT.println("correctionEnergy = " + correctionEnergy);
        val correctionEnergy = 0.0;
        //Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val totalEnergy = directEnergy + reciprocalEnergy + (correctionEnergy + selfEnergy);

        timer.stop(TIMER_INDEX_TOTAL);
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
        finish ateach(p in atoms) {
            val localAtoms = atoms(p);
            for (l in 0..(localAtoms.size-1)) {
                val atom = localAtoms(l);
                val centre = atom.centre;
                val charge = atom.charge;
                // get subcell i,j,k
                val i = (centre.i / halfCutoff) as Int;
                val j = (centre.j / halfCutoff) as Int;
                val k = (centre.k / halfCutoff) as Int;
                async at(subCellsTemp.dist(i,j,k)) {
                    atomic subCellsTemp(i,j,k).add(new PointCharge(centre, charge));
                }
            }
        }
        val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach([i,j,k] in subCells) {
            subCells(i,j,k) = subCellsTemp(i,j,k).toArray();
        }
    }

    /**
     * At each place, fetch all required atoms from neighbouring
     * places for direct calculation.
     */
    private def prefetchRemoteAtoms() {
        timer.start(TIMER_INDEX_PREFETCH);
		val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
		val atomsCache = this.atomsCache; // TODO shouldn't be necessary XTENLANG-1913
        finish ateach(p in atomsCache) {
            val myAtomsCache = atomsCache(here.id);
            if (myAtomsCache != null) {
                val haloPlaces = new HashMap[Int,ArrayList[Point(3)]](8); // a place may have up to 8 immediate neighbours in the two block-divided dimensions
                
                // separate the halo subcells into partial lists stored at each nearby place
                for (boxIndex in myAtomsCache.region) {
                    val placeId = subCells.dist(boxIndex).id;
                    if (placeId != here.id) {
                        var haloForPlace : ArrayList[Point(3)] = haloPlaces.getOrElse(placeId, null);
                        if (haloForPlace == null) {
                            haloForPlace = new ArrayList[Point(3)]();
                            haloPlaces.put(placeId, haloForPlace);
                        }
                        haloForPlace.add(boxIndex);
                    }
                }

                // retrieve the partial list for each place and store into my LET
                finish for (placeEntry in haloPlaces.entries()) async {
                    val placeId = placeEntry.getKey();
                    val haloForPlace = placeEntry.getValue();
                    val haloListArray = haloForPlace.toArray();
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
        val atoms = new Array[Rail[PointCharge]](boxList.size, 
                                                 (i : Int) => subCells(boxList(i)));
        return atoms;
    }

    public def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        val cutoffSquared = cutoff*cutoff;
		val numSubCells = this.numSubCells; // TODO shouldn't be necessary XTENLANG-1913
		val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
		val atomsCache = this.atomsCache; // TODO shouldn't be necessary XTENLANG-1913
		val imageTranslations = this.imageTranslations; // TODO shouldn't be necessary XTENLANG-1913
		val beta = this.beta; // TODO shouldn't be necessary XTENLANG-1913
        val directEnergy = finish(SumReducer()) {
            ateach(place in Dist.makeUnique()) {
                val remoteAtoms = atomsCache(here.id);
                val translations = imageTranslations(here.id);
                val localSubCells = subCells.getLocalPortion();
                val localRegion = localSubCells.region as Region(3){rect};
                for (p in localRegion) async {
                    val thisCell = localSubCells(p) as Rail[PointCharge];
                    var myDirectEnergy : Double = 0.0;
                    for (var i : Int = p(0)-2; i<=p(0); i++) {
                        var n1 : Int = 0;
                        if (i < 0) {
                            n1 = -1;
                        } // can't have (i > numSubCells+1)
                        for (var j : Int = p(1)-2; j<=p(1)+2; j++) {
                            var n2 : Int = 0;
                            if (j < 0) {
                                n2 = -1;
                            } else if (j > numSubCells-1) {
                                n2 = 1;
                            }
                            for (var k : Int = p(2)-2; k<=p(2)+2; k++) {
                                var n3 : Int = 0;
                                if (k < 0) {
                                    n3 = -1;
                                } else if (k > numSubCells-1) {
                                    n3 = 1;
                                }
                                // interact with "left half" of other boxes i.e. only boxes with i<=p(0)
                                if (i < p(0) || (i == p(0) && j < p(1)) || (i == p(0) && j == p(1) && k < p(2))) {
                                    val translation = translations(n1,n2,n3);
                                    val otherSubCellLocation = subCells.dist(i,j,k);
                                    val otherCell : Rail[PointCharge];
                                    if (otherSubCellLocation == here) {
                                        otherCell = subCells(i,j,k);
                                    } else {
                                        // other subcell is remote; use cached atoms
                                        otherCell = remoteAtoms(i,j,k);
                                    }
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
                                                myDirectEnergy += imageDirectComponent;
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
                                myDirectEnergy += directComponent;
                            }
                        }
                    }
                    offer myDirectEnergy;
                }
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
		val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
        val selfEnergy = finish(SumReducer()) {
            ateach(place in Dist.makeUnique()) {
                val localSubCells = subCells.getLocalPortion();
                val localRegion = localSubCells.region as Region(3){rect};
                var mySelfEnergy : Double = 0.0;
                for ([i,j,k] in localRegion) {
                    val thisCell = localSubCells(i,j,k) as Rail[PointCharge];
                    for (thisAtom in 0..(thisCell.size-1)) {
                        mySelfEnergy += thisCell(thisAtom).charge * thisCell(thisAtom).charge;
                    }
                }
                offer mySelfEnergy;
            }
        };
        timer.stop(TIMER_INDEX_SELF);
        return selfEnergy * -beta / Math.sqrt(Math.PI);
    }

    /** 
     * Calculates the gridded charge array Q as defined in Eq. 4.6,
     * using Cardinal B-spline interpolation.
     * Spline interpolation code inspired by GROMACS http://www.gromacs.org/
     */
    public def gridCharges() {
        timer.start(TIMER_INDEX_GRIDCHARGES);

        val gridSize = this.gridSize; // TODO shouldn't be necessary XTENLANG-1913
		val numSubCells = this.numSubCells; // TODO shouldn't be necessary XTENLANG-1913
		val splineOrder = this.splineOrder; // TODO shouldn't be necessary XTENLANG-1913
		val gridDist = this.gridDist; // TODO shouldn't be necessary XTENLANG-1913
		val subCells = this.subCells; // TODO shouldn't be necessary XTENLANG-1913
        val scalingVector = this.scalingVector; // TODO shouldn't be necessary XTENLANG-1913;
		val Q = this.Q; // TODO shouldn't be necessary XTENLANG-1913;
        finish ateach(place1 in Dist.makeUnique()) {
            val qLocal = Q.getLocalPortion() as Array[Complex](3){rect};
            qLocal.clear();
            val gridSize0 = gridSize(0);
            val gridSize2 = gridSize(2);
			val gridRegion = gridDist.region;
            val localSubCells = subCells.getLocalPortion();
            val place1Region = localSubCells.region as Region(3){rect};
            if (!place1Region.isEmpty()) {
                val place1HaloRegion = getGridRegionForSubcellAtoms(gridSize, numSubCells, splineOrder, place1Region);
                val myQ = new Array[Double](place1HaloRegion);
                val splines = new Array[Double](0..2 * 0..(splineOrder-1));

                for ([x,y,z] in place1Region) {
                    val thisCell = localSubCells(x,y,z) as Rail[PointCharge];
                    val atomOffset = new Array[Double](3);
                    for (atomIndex in 0..(thisCell.size-1)) {
                        val atom = thisCell(atomIndex);
                        val q = atom.charge;
                        val u = atom.centre.scale(scalingVector); // TODO general non-cubic
                        val u1c = Math.ceil(u.i) as Int;
                        val u2c = Math.ceil(u.j) as Int;
                        val u3c = Math.ceil(u.k) as Int;
                        atomOffset(0) = u1c - u.i;
                        atomOffset(1) = u2c - u.j;
                        atomOffset(2) = u3c - u.k;

                        for(j in 0..2) {
                            val offset = atomOffset(j);
                        
                            splines(j,splineOrder-1) = 0.0;
                            splines(j,1) = offset;
                            splines(j,0) = 1.0 - offset;
                        
                            for(k in 3..(splineOrder-1)) {
                                val div = 1.0 / (k-1.0);    
                                splines(j,k-1) = div * offset * splines(j,k-2);
                                for(l in 1..(k-2)) {
                                    splines(j,k-l-1) = div * ((offset+l) * splines(j,k-l-2) + (k-l-offset) * splines(j,k-l-1));
                                }
                                splines(j,0) = div * (1.0-offset) * splines(j,0);
                            }
                        
                            val div = 1.0 / (splineOrder-1);
                            splines(j,splineOrder-1) = div * offset * splines(j,splineOrder-2);
                            for(l in 1..(splineOrder-2)) {
                                splines(j,splineOrder-l-1) = div * 
                                    ((offset+l) * splines(j,splineOrder-l-2) + (splineOrder-l-offset) * splines(j,splineOrder-l-1));
                            }
                            splines(j,0) = div * (1.0-offset) * splines(j,0);
                        } 

                        for (i in 0..(splineOrder-1)) {
                            val k1 = (u1c - i - 1);
                            val iVal = q * splines(0, i);
                            for (j in 0..(splineOrder-1)) {
                                val k2 = (u2c - j - 1);
                                val jVal = iVal * splines(1, j);
                                for (k in 0..(splineOrder-1)) {
                                    val k3 = (u3c - k - 1);
                                    val kVal = jVal * splines(2, k);
                                    // because array is not divided in the z (k3) dimension, we can apply periodicity in that dimension 
                                    val wrapk3 = k3 < 0 ? (k3 + gridSize2) : (k3 >= gridSize2 ? (k3 - gridSize2) : k3);
                                    myQ(k1,k2,wrapk3) = myQ(k1,k2,wrapk3) + kVal;
                                }
                            }
                        }
                    }
                }

                // scatter myQ and accumulate to distributed Q
                finish for (place2 in gridDist.places()) {
                    val place2Region = gridDist.get(place2) as Region(3){rect};
                    // each halo region could be scattered to other regions in up to four chunks,
                    // as the region is periodic and divided along two dimensions.
                    val shiftX = (place1HaloRegion.max(0) < place2Region.min(0) || place1HaloRegion.min(0) < gridRegion.min(0)) ? gridSize0 : ((place1HaloRegion.min(0) > place2Region.max(0) || place1HaloRegion.max(0) > gridRegion.max(0)) ? -gridSize0 : 0);
                    val shiftY = (place1HaloRegion.max(1) < place2Region.min(1) || place1HaloRegion.min(1) < gridRegion.min(1)) ? gridSize0 : ((place1HaloRegion.min(1) > place2Region.max(1) || place1HaloRegion.max(1) > gridRegion.max(1)) ? -gridSize0 : 0);
                    scatterAndReduceShiftedGridContribution(myQ, Point.make(0,0,0), place2Region, Q, place2);
                    if (shiftX != 0) {
                        scatterAndReduceShiftedGridContribution(myQ, Point.make(shiftX,0,0), place2Region, Q, place2);
                    }
                    if (shiftY != 0) {
                        scatterAndReduceShiftedGridContribution(myQ, Point.make(0,shiftY,0), place2Region, Q, place2);
                        if (shiftX != 0) {
                           scatterAndReduceShiftedGridContribution(myQ, Point.make(shiftX,shiftY,0), place2Region, Q, place2);  
                        }
                    }
                  
                }
            }
        }

        timer.stop(TIMER_INDEX_GRIDCHARGES);
    }


    /**
     * Scatters and reduces that portion of the sourceGrid array
     * that overlaps with the Q grid array at the target place,
     * when shifted by an amount.  The shift is required to represent
     * the <= 4 possible subregions that may be transferred from a
     * place to a "neighbouring" place, given that in the case of a 
     * periodic array divided in two dimensions between four places,
     * one place may "neighbour" another across corner diagonals in 
     * four different directions.
     * N.B. the case of a periodic array divided in three dimensions
     * is even worse!
     * TODO a more general solution to this problem with be required to solve
     * XTENLANG-1373 (in conjuction with XTENLANG-1365)
	 * TODO should use instance fields instead of all those parameters - XTENLANG-1913
     */
    private static def scatterAndReduceShiftedGridContribution(sourceGrid : Array[Double]{self.rect,self.rank==3},
                                             shift : Point(3),
										     targetRegion : Region(3){rect},
											 Q : DistArray[Complex](3),
                                             targetPlace : Place) {
        val overlapRegion = (sourceGrid.region + shift && targetRegion) as Region(3){rect};
        if (! overlapRegion.isEmpty()) {
            val sourceRegion = (overlapRegion - shift) as Region(3){rect};
            val overlap = new Array[Double](overlapRegion.size());

            var l : Int = 0;
            for ([i,j,k] in sourceRegion) {
                overlap(l++) = sourceGrid(i,j,k);
            }
            async at(targetPlace) {
                val localQ = Q.getLocalPortion();
                atomic {
                    var m : Int = 0;
                    for ([i,j,k] in overlapRegion) {
                        localQ(i,j,k) = localQ(i,j,k) + overlap(m++);
                    }
                }
            }
        }
    }

    /**
     * Given the minimum and maximum boundaries of a (rectangular) region
     * of subcells, returns the subregion which is contributed to by the atoms
     * in those subcells.
     * Each subcell only contributes to a portion of the grid, depending on
     * the spline order.  (A larger spline order spreads the atoms in a subcell
     * over a larger area of the grid.)
	 * TODO general non-cubic grid
	 * TODO should use instance fields instead of all those parameters - XTENLANG-1913
     */
    private static @Inline def getGridRegionForSubcellAtoms(gridSize : Rail[Int], numSubCells : Int, splineOrder : Int, subCellRegion : Region(3){rect}) {
        val gridPointsPerSubCell = (gridSize(0) as Double) / numSubCells;
        val minX = Math.floor(subCellRegion.min(0) * gridPointsPerSubCell) as Int - splineOrder + 1;
        val maxX = Math.ceil((subCellRegion.max(0)+1) * gridPointsPerSubCell) as Int - 1;
        val minY = Math.floor(subCellRegion.min(1) * gridPointsPerSubCell) as Int - splineOrder + 1;
        val maxY = Math.ceil((subCellRegion.max(1)+1) * gridPointsPerSubCell) as Int - 1;
        return (minX..maxX) * (minY..maxY) * (0..(gridSize(2)-1));
    }

    /**
     * @return the approximation to the reciprocal energy ~E_rec as defined in Eq. 4.7
     */
    private def getReciprocalEnergy(thetaRecConvQ : DistArray[Complex]{self.dist==gridDist}) {
        timer.start(TIMER_INDEX_RECIPROCAL);

        val reciprocalEnergy = finish(SumReducer()) {
			val gridDist = this.gridDist; // TODO shouldn't be necessary XTENLANG-1913
			val Q = this.Q; // TODO shouldn't be necessary XTENLANG-1913
			//val thetaRecConvQ = this.thetaRecConvQ; // TODO shouldn't be necessary XTENLANG-1913
            ateach(place in Dist.makeUnique()) {
                var myReciprocalEnergy : Double = 0.0;
                val localQ = Q.getLocalPortion();
                val localThetaRecConvQ = thetaRecConvQ.getLocalPortion();
                val localRegion = localQ.region as Region(3){rect};
                for ([i,j,k] in localRegion) {
                    val gridPointContribution = localQ(i,j,k) * localThetaRecConvQ(i,j,k);
                    myReciprocalEnergy += gridPointContribution.re;
                }
                offer myReciprocalEnergy;
            }
        };

        timer.stop(TIMER_INDEX_RECIPROCAL);
        return reciprocalEnergy / 2.0;
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
        if (n == 4) {
            return bSpline4(u);
        } else if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 2) {
            return 1.0 - Math.abs(u - 1.0);
        } else {
            return u / (n - 1) * bSpline(n-1, u) + (n - u) / (n - 1) * bSpline(n-1, u-1.0);
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
    
    /** Gets scaled fractional coordinate u as per Eq. 3.1 - general rectangular */
    public static @Inline def getScaledFractionalCoordinates(K1 : Double, K2 : Double, K3 : Double, edgeReciprocals : Rail[Vector3d], r : Vector3d) : Vector3d {
        // this method allows general non-rectangular cells
        return Vector3d(edgeReciprocals(0).mul(K1).dot(r), edgeReciprocals(1).mul(K2).dot(r), edgeReciprocals(2).mul(K3).dot(r));
    }

    /**
     * Gets the volume V of the unit cell.
     */
    private def getVolume() {
        return edges(0).cross(edges(1)).dot(edges(2));
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
