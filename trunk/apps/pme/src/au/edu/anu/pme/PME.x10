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
package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.fft.Distributed3dFft;
import au.edu.anu.util.Timer;
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
    private val gridSize : Array[Int]{rail};

    /** Double representations of the various grid dimensions */
    private val K1 : Double;
    private val K2 : Double; 
    private val K3 : Double;

    /** The edges of the unit cell. */
    private val edges : Array[Vector3d]{rail};
    private val edgeLengths : Array[Double]{rail};

    /** The conjugate reciprocal vectors for each dimension. */
    private val edgeReciprocals : Array[Vector3d]{rail};

    private val gridRegion : Region{rank==3,rect};
    private val gridDist : Dist{rank==3,rect};
    
    /** The order of B-spline interpolation */
    private val splineOrder : Int;

    /** The Ewald coefficient beta */
    private val beta : Double;

    /** The direct sum cutoff distance in Angstroms */
    private val cutoff : Double;

    /** 
     * Translation vectors for neighbouring unit cells 
     * (the 26 cells surrounding the origin cell)
     * Dimensions are:
     * 0: a "virtual" dimension used to replicate the data across all places
     *   // TODO should be done as a global array XTENLANG-787
     * 1: x translation (difference between x-coordinate of sub-cells
     * 2: y translation
     * 3: z translation
     */
    private val imageTranslations : DistArray[Vector3d]{rank==4,rect};

    /** The atoms in the simulation, divided up into a distributed array of Arrays, one for each place. */
    private val atoms : DistArray[Array[MMAtom]{rail}]{rail};

    private val B : DistArray[Double]{self.dist==gridDist};
    private val C : DistArray[Double]{self.dist==gridDist};
    private val BdotC : DistArray[Double]{self.dist==gridDist};

    /** The reciprocal pair potential as defined in eq. 4.7 */
    private val thetaRecConvQ : DistArray[Complex]{self.dist==gridDist};

    /** The gridded charge array Q as defined in Eq. 4.6 */
    private val Q : DistArray[Complex]{self.dist==gridDist};

    /** The inverse DFT of the Q array.  TODO this should be a scoped local variable in getEnergy() XTENLANG-??? */
    private val Qinv : DistArray[Complex]{self.dist==gridDist};

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
    private val subCells : PeriodicDistArray[Array[MMAtom]{rail}]{rank==3,rect};
    /** The number of sub-cells per side of the unit cell. */
    private val numSubCells : Int;

    /** 
     * A cache of "packed" representations of atoms from subcells
     * stored at other places.  This is used to prefetch atom data
     * for direct energy calculation.
     */
    private val packedAtomsCache : DistArray[Array[Array[MMAtom.PackedRepresentation]{rail}]{rank==3,rect}]{rail};

    /**
     * Creates a new particle mesh Ewald method.
     * @param edges the edge vectors of the unit cell
     * @param gridSize the number of grid lines in each dimension of the unit cell
     * @param atoms the atoms in the unit cell
     * @param splineOrder the order n of B-spline interpolationb
     * @param beta the Ewald coefficient beta
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(edges : Array[Vector3d]{rail},
            gridSize : Array[Int]{rail},
            atoms: DistArray[Array[MMAtom]{rail}]{rail},
            splineOrder : Int,
            beta : Double,
            cutoff : Double) {
        this.gridSize = gridSize;
        K1 = gridSize(0) as Double;
        K2 = gridSize(1) as Double;
        K3 = gridSize(2) as Double;
        this.edges = edges;
        this.edgeLengths = new Array[Double](3, (i : Int) => edges(i).length());
        this.edgeReciprocals = new Array[Vector3d](3, (i : Int) => edges(i).inverse());

        this.atoms = atoms;
        gridRegion = (0..gridSize(0)-1) * (0..gridSize(1)-1) * (0..gridSize(2)-1);
        gridDist = Dist.makeBlockBlock(gridRegion, 0, 1);
        this.splineOrder = splineOrder;
        this.beta = beta;
        this.cutoff = cutoff;
        val imageTranslationRegion = (0..Place.MAX_PLACES-1) * (-1..1) * (-1..1) * (-1..1);
        this.imageTranslations = DistArray.make[Vector3d](Dist.makeBlock(imageTranslationRegion, 0), ([place,i,j,k] : Point(4)) => (edges(0).mul(i)).add(edges(1).mul(j)).add(edges(2).mul(k)));

        if (edgeLengths(0) % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + edgeLengths(0) + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(edgeLengths(0) / (cutoff/2.0)) as Int;
        val subCellRegion = (0..numSubCells-1) * (0..numSubCells-1) * (0..numSubCells-1);
        val subCells = PeriodicDistArray.make[Array[MMAtom]{rail}](Dist.makeBlockBlock(subCellRegion, 0, 1));
        Console.OUT.println("subCells dist = " + subCells.dist);
        this.subCells = subCells;
        this.numSubCells = numSubCells;

        packedAtomsCache = DistArray.make[Array[Array[MMAtom.PackedRepresentation]{rail}]{rank==3,rect}](Dist.makeUnique());
        finish ateach (p in packedAtomsCache) {
            val mySubCellRegion = (subCells.dist | here).region;
            if (! mySubCellRegion.isEmpty()) {
                val directRequiredRegion = ((mySubCellRegion.min(0) - 2)..(mySubCellRegion.max(0) + 2))
                                         * ((mySubCellRegion.min(1) - 2)..(mySubCellRegion.max(1) + 2))
                                         * ((mySubCellRegion.min(2) - 2)..(mySubCellRegion.max(2) + 2));
                packedAtomsCache(p) = new Array[Array[MMAtom.PackedRepresentation]{rail}](directRequiredRegion);
            }
        }


        Console.OUT.println("gridDist = " + gridDist);

        // TODO should be able to automatically zero arrays of Complex
        Q = DistArray.make[Complex](gridDist, (Point) => Complex.ZERO);
        BdotC = DistArray.make[Double](gridDist);
        thetaRecConvQ = DistArray.make[Complex](gridDist);

        // TODO following arrays should be scoped local variables XTENLANG-???
        Qinv = DistArray.make[Complex](gridDist, (Point) => Complex.ZERO);
        temp = DistArray.make[Complex](gridDist, (Point) => Complex.ZERO);
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
        finish for (place1 in gridDist.places()) async at(place1) {
            for (p in gridDist.get(here)) {
                BdotC(p) = B(p) * C(p);
            }
        }
        divideAtomsIntoSubCells();
        timer.stop(TIMER_INDEX_SETUP);
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        finish {
            async { prefetchPackedAtoms(); }

            gridCharges();


            timer.start(TIMER_INDEX_INVFFT);
            new Distributed3dFft(gridSize(0), Q, Qinv, temp).doFFT3d(false);
            timer.stop(TIMER_INDEX_INVFFT);

            timer.start(TIMER_INDEX_THETARECCONVQ);
            // create F^-1(thetaRecConvQ)
            finish ateach (place in Dist.makeUnique()) {
                val myGridRegion = gridDist.get(here) as Region{rank==3,rect};
                for ([i,j,k] in myGridRegion) {
                    thetaRecConvQ(i,j,k) = BdotC(i,j,k) * Qinv(i,j,k);
                }
            }
            // and do inverse FFT
            new Distributed3dFft(gridSize(0), thetaRecConvQ, thetaRecConvQ, temp).doFFT3d(true);
            timer.stop(TIMER_INDEX_THETARECCONVQ);
        }

        val reciprocalEnergy = getReciprocalEnergy();
        val selfEnergy = getSelfEnergy();
        val directEnergy = getDirectEnergy();

        Console.OUT.println("directEnergy = " + directEnergy);
        Console.OUT.println("selfEnergy = " + selfEnergy);
        //Console.OUT.println("correctionEnergy = " + correctionEnergy);
        val correctionEnergy = 0.0;
        Console.OUT.println("reciprocalEnergy = " + reciprocalEnergy);
        val totalEnergy = directEnergy + reciprocalEnergy + (correctionEnergy + selfEnergy);

        timer.stop(TIMER_INDEX_TOTAL);
        return totalEnergy;
    }

    private def divideAtomsIntoSubCells() {
        val subCellsTemp = DistArray.make[ArrayList[MMAtom]](subCells.dist, (Point) => new ArrayList[MMAtom]());
        finish ateach (p1 in atoms) {
            val localAtoms = atoms(p1);
            finish for ([i] in 0..localAtoms.size()-1) async {
                val atom = localAtoms(i);
                val charge = atom.charge;
                val centre = atom.centre;
                val index = getSubCellIndex(atom);
                at (subCellsTemp.dist(index)) {
                    val remoteAtom = new MMAtom(centre, charge);
                    atomic{subCellsTemp(index).add(remoteAtom);}
                }
            }
        }
        finish ateach (p in subCells.dist) {
            subCells(p) = subCellsTemp(p).toArray();
        }
    }

    /**
     * Gets the index of the sub-cell for this atom within
     * the grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def getSubCellIndex(atom : MMAtom) : Point(3) {
        // TODO assumes cubic unit cell
        val halfCutoff = (cutoff / 2.0);
        val r = Vector3d(atom.centre);
        val i = (r.i / halfCutoff) as Int;
        val j = (r.j / halfCutoff) as Int;
        val k = (r.k / halfCutoff) as Int;
        return Point.make(i, j, k);
    }

    /**
     * At each place, fetch all required atoms from neighbouring
     * places for direct calculation.
     */
    private def prefetchPackedAtoms() {
        timer.start(TIMER_INDEX_PREFETCH);
        finish for (place in subCells.dist.places()) async at (place) {
            val myPackedAtoms = packedAtomsCache(here.id);
            val haloPlaces = new HashMap[Int,ArrayList[Point(3)]](8); // a place may have up to 8 immediate neighbours in the two block-divided dimensions
            
            // separate the halo subcells into partial lists stored at each nearby place
            for (boxIndex in myPackedAtoms.region) {
                val placeId = subCells.periodicDist(boxIndex).id;
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
                val packedForPlace = at (Place.place(placeId)) { getPackedAtomsForSubcellList(haloListArray)};
                for ([i] in 0..haloListArray.size-1) {
                    myPackedAtoms(haloListArray(i)) = packedForPlace(i);
                }
            }
        }
        timer.stop(TIMER_INDEX_PREFETCH);
    }

    /**
     * Given a list of subcell indices as Point(3) stored at a single
     * place, returns an Array, each element of which is in turn
     * a Array of MMAtom.PackedRepresentation containing the 
     * packed atoms for each subcell.
     */
    private def getPackedAtomsForSubcellList(boxList : Array[Point(3)]{rail}) {
        val packedAtomList = new Array[Array[MMAtom.PackedRepresentation]{rail}](boxList.size, 
                                                            (i : Int) => 
                                                                getPackedAtomsForSubCell(boxList(i))
                                                            );
        return packedAtomList;
    }

    public def getDirectEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        val cutoffSquared = cutoff*cutoff;
        val directEnergy = finish(SumReducer()) {
            ateach (p in subCells) {
                var myDirectEnergy : Double = 0.0;
                val thisCell = subCells(p);
                val packedAtoms = packedAtomsCache(here.id);
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
                                val translation = imageTranslations(here.id,n1,n2,n3);
                                val otherSubCellLocation = subCells.periodicDist(i,j,k);
                                if (otherSubCellLocation == here) {
                                    val otherCell = subCells(i,j,k);
                                    for ([otherAtom] in 0..otherCell.size-1) {
                                        val imageLoc = otherCell(otherAtom).centre + translation;
                                        for ([thisAtom] in 0..thisCell.size-1) {
                                            val rSquared = thisCell(thisAtom).centre.distanceSquared(imageLoc);
                                            if (rSquared < cutoffSquared) {
                                                val r = Math.sqrt(rSquared);
                                                val chargeProduct = thisCell(thisAtom).charge * otherCell(otherAtom).charge;
                                                val imageDirectComponent = chargeProduct * Math.erfc(beta * r) / r;
                                                myDirectEnergy += imageDirectComponent;
                                            }
                                        }
                                    }
                                } else {
                                    // other subcell is remote; use cached packed atoms
                                    val otherCellPacked = packedAtoms(i,j,k);
                                    for ([otherAtom] in 0..otherCellPacked.size-1) {
                                        val imageLoc = otherCellPacked(otherAtom).centre + translation;
                                        for ([thisAtom] in 0..thisCell.size-1) {
                                            val rSquared = thisCell(thisAtom).centre.distanceSquared(imageLoc);
                                            if (rSquared < cutoffSquared) {
                                                val r = Math.sqrt(rSquared);
                                                val chargeProduct = thisCell(thisAtom).charge * otherCellPacked(otherAtom).charge;
                                                val imageDirectComponent = chargeProduct * Math.erfc(beta * r) / r;
                                                myDirectEnergy += imageDirectComponent;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // atoms in same cell
                for ([i] in 0..thisCell.size-1) {
                    for ([j] in 0..i-1) {
                        val rjri = thisCell(j).centre - thisCell(i).centre;
                        val rSquared = rjri.lengthSquared();
                        if (rSquared < cutoffSquared) {
                            val r = Math.sqrt(rSquared);
                            val directComponent = thisCell(i).charge * thisCell(j).charge * Math.erfc(beta * r) / r;
                            myDirectEnergy += directComponent;
                        }
                    }
                }
                offer myDirectEnergy;
            }
        };
        
        timer.stop(TIMER_INDEX_DIRECT);
        return directEnergy;
    }

    /*
     * Returns atom charges and coordinates for a sub-cell, in packed representation
     */
    private def getPackedAtomsForSubCell(subCellIndex : Point(3)) : Array[MMAtom.PackedRepresentation]{rail} {
        val subCell = subCells(subCellIndex);
        return new Array[MMAtom.PackedRepresentation](subCell.size, (i : Int) => subCell(i).getPackedRepresentation());
    }

    /**
     * Returns the self energy as defined in Eq. 2.5, which is the
     * contribution of the interaction in reciprocal space of each
     * atom with itself.  This is subtracted from the final energy.
     */
    public def getSelfEnergy() : Double {
        timer.start(TIMER_INDEX_SELF);
        val selfEnergy = finish(SumReducer()) {
            ateach (p in subCells.dist) {
                val thisCell = subCells(p);
                var mySelfEnergy : Double = 0.0;
                for ([thisAtom] in 0..thisCell.size-1) {
                    mySelfEnergy += thisCell(thisAtom).charge * thisCell(thisAtom).charge;
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
     */
    public def gridCharges() {
        timer.start(TIMER_INDEX_GRIDCHARGES);
        finish ateach (place1 in Dist.makeUnique()) {
            val place1Region = subCells.dist.get(here);
            if (!place1Region.isEmpty()) {
                val place1HaloRegion = getGridRegionForSubcellAtoms(place1Region);
                val myQ = new Array[Double](place1HaloRegion, (Point) => 0.0);

                for (p in place1Region) {
                    val thisCell = subCells(p);
                    for (atomIndex in thisCell) {
                        val atom = thisCell(atomIndex);
                        val q = atom.charge;
                        val u = getScaledFractionalCoordinates(atom.centre); // TODO general non-cubic
                        val u1c = Math.ceil(u.i) as Int;
                        val u2c = Math.ceil(u.j) as Int;
                        val u3c = Math.ceil(u.k) as Int;
                        for ([i] in 1..splineOrder) {
                            val k1 = (u1c - i);
                            for ([j] in 1..splineOrder) {
                                val k2 = (u2c - j);
                                for ([k] in 1..splineOrder) {
                                    val k3 = (u3c - k);
                                    val gridPointContribution = q
                                                 * bSpline(splineOrder, u.i - k1)
                                                 * bSpline(splineOrder, u.j - k2)
                                                 * bSpline(splineOrder, u.k - k3);
                                    // because array is not divided in the z (k3) dimension, we can apply periodicity in that dimension 
                                    val wrapk3 = k3 < 0 ? (k3 + gridSize(2)) : (k3 >= gridSize(2) ? (k3 - gridSize(2)) : k3);
                                    myQ(k1,k2,wrapk3) = myQ(k1,k2,wrapk3) + gridPointContribution;
                                }
                            }
                        }
                    }
                }

                // scatter myQ and accumulate to distributed Q
                finish for (place2 in gridDist.places()) {
                    val place2Region = gridDist.get(place2);
                    // each halo region could be scattered to other regions in up to four chunks,
                    // as the region is periodic and divided along two dimensions.
                    val shiftX = (place1HaloRegion.max(0) < place2Region.min(0) || place1HaloRegion.min(0) < gridRegion.min(0)) ? gridSize(0) : ((place1HaloRegion.min(0) > place2Region.max(0) || place1HaloRegion.max(0) > gridRegion.max(0)) ? -gridSize(0) : 0);
                    val shiftY = (place1HaloRegion.max(1) < place2Region.min(1) || place1HaloRegion.min(1) < gridRegion.min(1)) ? gridSize(0) : ((place1HaloRegion.min(1) > place2Region.max(1) || place1HaloRegion.max(1) > gridRegion.max(1)) ? -gridSize(0) : 0);
                    scatterAndReduceShiftedGridContribution(myQ, Point.make(0,0,0), place2);
                    if (shiftX != 0) {
                        scatterAndReduceShiftedGridContribution(myQ, Point.make(shiftX,0,0), place2);
                    }
                    if (shiftY != 0) {
                        scatterAndReduceShiftedGridContribution(myQ, Point.make(0,shiftY,0), place2);
                        if (shiftX != 0) {
                           scatterAndReduceShiftedGridContribution(myQ, Point.make(shiftX,shiftY,0), place2);  
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
     */
    private def scatterAndReduceShiftedGridContribution(sourceGrid : Array[Double]{self.rect,self.rank==3},
                                             shift : Point(3),
                                             targetPlace : Place) {
        val targetRegion = (gridDist | targetPlace).region;
        val overlapRegion = (sourceGrid.region + shift && targetRegion) as Region{rank==3,rect};
        val sourceRegion = (overlapRegion - shift) as Region{rank==3,rect};
        if (! overlapRegion.isEmpty()) {
            val overlap = new Array[Double](overlapRegion.size());

            var i : Int = 0;
            for (p in sourceRegion) {
                overlap(i++) = sourceGrid(p);
            }
            at (targetPlace) {
                atomic {
                    var j : Int = 0;
                    for (p in overlapRegion) {
                        Q(p) = Q(p) + overlap(j++);
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
     */
    private def getGridRegionForSubcellAtoms(subCellRegion : Region{rank==3}) {
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
    private def getReciprocalEnergy() {
        timer.start(TIMER_INDEX_RECIPROCAL);

        val reciprocalEnergy = finish(SumReducer()) {
            for (place1 in gridDist.places()) async at(place1) {
                var myReciprocalEnergy : Double = 0.0;
                for (p in gridDist | here) {
                    val gridPointContribution = Q(p) * thetaRecConvQ(p);
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
        finish ateach (place1 in Dist.makeUnique()) {
            for ([m1,m2,m3] in B.dist(here)) {
                val m1D = m1 as Double;
                var sumK1 : Complex = Complex.ZERO;
                for ([k] in 0..(splineOrder-2)) {
                    sumK1 = sumK1 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m1D * k / K1 * Complex.I);
                }
                val b1 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m1D / K1 * Complex.I) / sumK1).abs();

                val m2D = m2 as Double;
                var sumK2 : Complex = Complex.ZERO;
                for ([k] in 0..(splineOrder-2)) {
                    sumK2 = sumK2 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m2D * k / K2 * Complex.I);
                }
                val b2 = (Math.exp(2.0 * Math.PI * (splineOrder-1) * m2D / K2 * Complex.I) / sumK2).abs();
                    
                val m3D = m3 as Double;
                var sumK3 : Complex = Complex.ZERO;
                for ([k] in 0..(splineOrder-2)) {
                    sumK3 = sumK3 + bSpline(splineOrder, k+1) * Math.exp(2.0 * Math.PI * m3D * k / K3 * Complex.I);
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
        finish ateach (place1 in Dist.makeUnique()) {
            for ([m1,m2,m3] in C.dist(here)) {
                val m1prime = m1 <= K1/2 ? m1 : m1 - K1;
                val m2prime = m2 <= K2/2 ? m2 : m2 - K2;
                val m3prime = m3 <= K3/2 ? m3 : m3 - K3;
                val mVec = edgeReciprocals(0).mul(m1prime).add(edgeReciprocals(1).mul(m2prime)).add(edgeReciprocals(2).mul(m3prime));
                val mSquared = mVec.dot(mVec);
                C(m1,m2,m3) = Math.exp(-(Math.PI*Math.PI) * mSquared / (beta * beta)) / (mSquared * Math.PI * V);
            }
        }
        at (C.dist(0,0,0)) {
            C(0,0,0) = 0.0;
        }
    }

    /* 
     * Gets the nth order B-spline M_n(u) as per Eq. 4.1
     */
    public final static def bSpline(n : Int, u : Double) : Double {
        if (u < 0.0 || u > n) {
            return 0.0;
        } else if (n == 4) {
            return bSpline4(u);
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

    /** Gets scaled fractional coordinate u as per Eq. 3.1 - cubic only */
    public def getScaledFractionalCoordinates(r : Point3d) : Vector3d {
        return Vector3d(edgeReciprocals(0).i * K1 * r.i, edgeReciprocals(1).j * K2 * r.j, edgeReciprocals(2).k * K3 * r.k);
    }
    
    /** Gets scaled fractional coordinate u as per Eq. 3.1 - general rectangular */
    public def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
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
        public def apply(a:Double, b:Double) = (a + b);
    }
}
