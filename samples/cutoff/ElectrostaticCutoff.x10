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
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.PointCharge;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.util.Timer;
import au.edu.anu.chem.mm.TestElectrostatic;

import x10.util.ArrayList;
import x10.util.HashMap;
import x10.util.Pair;

/**
 * This class implements electrostatic potential calculation with cutoff
 * and periodic boundary conditions.
 */
public class ElectrostaticCutoff extends TestElectrostatic {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL : Int = 0;
    public static val TIMER_INDEX_DIRECT : Int = 1;
    public static val TIMER_INDEX_PREFETCH : Int = 2;
    public static val TIMER_INDEX_SETUP : Int = 3;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(3);

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
     * @param edges the edge vectors of the unit cell
     * @param gridSize the number of grid lines in each dimension of the unit cell
     * @param atoms the atoms in the unit cell
     * @param splineOrder the order n of B-spline interpolationb
     * @param beta the Ewald coefficient beta
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(size:Double,
            cutoff : Double) {

        this.cutoff = cutoff;
        this.imageTranslations = DistArray.make[Array[Vector3d](3){rect}](
            Dist.makeUnique(), 
            new Array[Vector3d](-1..1 * -1..1 * -1..1, 
                ([i,j,k] : Point(3)) => Vector3d(i*size, j*size, k*size)) 
        );

        if (cutoff > size) {
            throw new IllegalArgumentException("error: cutoff " + cutoff + " is greater than size " + size);
        }
        if (size % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + size + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(size / (cutoff/2.0)) as Int;
        val subCellRegion = 0..(numSubCells-1) * 0..(numSubCells-1) * 0..(numSubCells-1);
        val subCells = DistArray.make[Rail[PointCharge]](new PeriodicDist(Dist.makeBlockBlock(subCellRegion, 0, 1)));
        //Console.OUT.println("subCells dist = " + subCells.dist);
        this.subCells = subCells;
        this.numSubCells = numSubCells;

        val atomsCache = DistArray.make[Array[Rail[PointCharge]]{rank==3,rect}](Dist.makeUnique());
        finish ateach(p in atomsCache) {
            val mySubCellRegion = subCells.dist(here);
            if (! mySubCellRegion.isEmpty()) {
                val directRequiredRegion = ((mySubCellRegion.min(0) - 2)..(mySubCellRegion.max(0) + 2))
                                         * ((mySubCellRegion.min(1) - 2)..(mySubCellRegion.max(1) + 2))
                                         * ((mySubCellRegion.min(2) - 2)..(mySubCellRegion.max(2) + 2));
                atomsCache(p) = new Array[Rail[PointCharge]](directRequiredRegion);
            }
        }
        this.atomsCache = atomsCache;
    }

    /**
     * Divide the atoms into a grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def assignAtomsToSubCells(atoms: DistArray[Rail[MMAtom]](1)) {
        timer.start(TIMER_INDEX_SETUP);
        val halfCutoff = (cutoff / 2.0);
        val subCellsTemp = DistArray.make[ArrayList[PointCharge]](subCells.dist, (Point) => new ArrayList[PointCharge]());
        val halfNumSubCells = this.numSubCells / 2;
        finish ateach(p in atoms) {
            val localAtoms = atoms(p);
            for (l in 0..(localAtoms.size-1)) {
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
            subCells(i,j,k) = subCellsTemp(i,j,k).toArray();
        }
        timer.stop(TIMER_INDEX_SETUP);
    }

    public def getEnergy() {
        prefetchRemoteAtoms();
        return getDirectEnergy();
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
                    val haloListArray = haloForPlace.toArray();
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
		val subCellsDist = this.subCells.dist;
		val numSubCells = this.numSubCells; // TODO shouldn't be necessary XTENLANG-1913
		val atomsCache = this.atomsCache; // TODO shouldn't be necessary XTENLANG-1913
		val imageTranslations = this.imageTranslations; // TODO shouldn't be necessary XTENLANG-1913
        val directEnergy = new Accumulator[Double](Reducible.SumReducer[Double]());
        finish ateach(place in Dist.makeUnique()) {
            val cachedAtoms = atomsCache(here.id);
            val translations = imageTranslations(here.id);
            val localRegion = subCellsDist(here) as Region(3){rect};
            for ([x,y,z] in localRegion) async {
                val thisCell = cachedAtoms(x,y,z) as Rail[PointCharge];
                var myDirectEnergy : Double = 0.0;
                for (var i : Int = x-2; i<=x; i++) {
                    var n1 : Int = 0;
                    if (i < 0) {
                        n1 = -1;
                    } // can't have (i > numSubCells+1)
                    for (var j : Int = y-2; j<=y+2; j++) {
                        var n2 : Int = 0;
                        if (j < 0) {
                            n2 = -1;
                        } else if (j > numSubCells-1) {
                            n2 = 1;
                        }
                        for (var k : Int = z-2; k<=z+2; k++) {
                            var n3 : Int = 0;
                            if (k < 0) {
                                n3 = -1;
                            } else if (k > numSubCells-1) {
                                n3 = 1;
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
                                            val imageDirectComponent = chargeProduct / r;
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
                            val directComponent = thisAtom.charge * otherAtom.charge / r;
                            myDirectEnergy += directComponent;
                        }
                    }
                }
                directEnergy <- myDirectEnergy;
            }
        }
        
        timer.stop(TIMER_INDEX_DIRECT);
        return directEnergy();
    }

    public static def main(args : Rail[String]) {
        val size:Double = 80.0;
        var numAtoms : Int;
        var cutoff : Double = 10.0;
        if (args.size > 0) {
            numAtoms = Int.parseInt(args(0));
            if (args.size > 1) {
                cutoff = Double.parseDouble(args(1));
            }
        } else {
            Console.ERR.println("usage: cutoff numAtoms [cutoff]");
            return;
        }
        Console.OUT.println("electrostatics for " + numAtoms + " atoms with cutoff " + cutoff + " (total size " + size + ")");
        new ElectrostaticCutoff(size, cutoff).test(numAtoms);
    }

    public def test(numAtoms:Int) {
        val atoms = generateAtoms(numAtoms);
        assignAtomsToSubCells(atoms);
        Console.OUT.println("Total PE: " + getEnergy());
        logTime("Setup",             ElectrostaticCutoff.TIMER_INDEX_SETUP,         timer);
        logTime("Prefetch",          ElectrostaticCutoff.TIMER_INDEX_PREFETCH,      timer);
        logTime("Direct",            ElectrostaticCutoff.TIMER_INDEX_DIRECT,        timer);
    }

}
