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
import au.edu.anu.util.Timer;

import x10.array.Array_3;
import x10.regionarray.Array;
import x10.regionarray.Region;
import x10.util.ArrayList;
import x10.util.Random;

/**
 * This class implements electrostatic potential calculation with cutoff
 * and periodic boundary conditions.
 */
public class ElectrostaticCutoff {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL:Long = 0;
    public static val TIMER_INDEX_DIRECT:Long = 1;
    public static val TIMER_INDEX_SETUP:Long = 2;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(3);

    /** The side length of the cubic simulation space. */
    private val size:Double;

    /** The direct sum cutoff distance in Angstroms */
    private val cutoff:Double;

    /** 
     * Translation vectors for neighbouring unit cells 
     * (the 26 cells surrounding the origin cell).
     * Dimensions are:
     * 0: x translation (difference between x-coordinate of sub-cells
     * 1: y translation
     * 2: z translation
     */
    private val imageTranslations:Array[Vector3d]{rect};


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
    private val subCells:Array_3[Rail[PointCharge]];

    /** The number of sub-cells per side of the unit cell. */
    private val numSubCells:Long;

    /**
     * @param size the length of a side of the cubic simulation space
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(size:Double, cutoff : Double) {
        this.size = size;
        this.cutoff = cutoff;
        this.imageTranslations = new Array[Vector3d](Region.makeRectangular(-1..1, -1..1, -1..1), 
                ([i,j,k] : Point(3)) => Vector3d(i*size, j*size, k*size));

        if (cutoff > size) {
            throw new IllegalArgumentException("error: cutoff " + cutoff + " is greater than size " + size);
        }
        if (size % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + size + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(size / (cutoff/2.0)) as Int;
        val subCells = new Array_3[Rail[PointCharge]](numSubCells, numSubCells, numSubCells);
        this.subCells = subCells;
        this.numSubCells = numSubCells;
    }

    /**
     * Divide the atoms into a grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def assignAtomsToSubCells(atoms:Rail[PointCharge]) {
        timer.start(TIMER_INDEX_SETUP);
        val halfCutoff = (cutoff / 2.0);
        val subCellsTemp = new Array_3[ArrayList[PointCharge]](subCells.numElems_1, subCells.numElems_2, subCells.numElems_3, (i:Long,j:Long,k:Long) => new ArrayList[PointCharge]());
        val halfNumSubCells = this.numSubCells / 2;
        for (l in 0..(atoms.size-1)) {
            val atom = atoms(l);
            val centre = atom.centre;
            // get subcell i,j,k
            val i = (centre.i / halfCutoff) as Long + halfNumSubCells;
            val j = (centre.j / halfCutoff) as Long + halfNumSubCells;
            val k = (centre.k / halfCutoff) as Long + halfNumSubCells;
            subCellsTemp(i,j,k).add(atom);
        }
        for([i,j,k] in subCells.indices()) {
            subCells(i,j,k) = subCellsTemp(i,j,k).toRail();
        }
        timer.stop(TIMER_INDEX_SETUP);
    }

    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        val cutoffSquared = cutoff*cutoff;
        val directEnergy = new Accumulator[Double](Reducible.SumReducer[Double]());
        finish for ([x,y,z] in subCells.indices()) async {
            val thisCell = subCells(x,y,z);
            var myDirectEnergy : Double = 0.0;
            for (i in (x-2)..x) {
                var n1:Long = 0;
                if (i < 0) {
                    n1 = -1;
                } // can't have (i > numSubCells+1)
                for (j in (y-2)..y) {
                    var n2:Long = 0;
                    if (j < 0) {
                        n2 = -1;
                    } else if (j > numSubCells-1) {
                        n2 = 1;
                    }
                    for (k in (z-2)..z) {
                        var n3:Long = 0;
                        if (k < 0) {
                            n3 = -1;
                        } else if (k > numSubCells-1) {
                            n3 = 1;
                        }
                        // interact with "left half" of other boxes i.e. only boxes with i<=x
                        if (i < x || (i == x && j < y) || (i == x && j == y && k < z)) {
                            val otherCell = subCells(i-n1*numSubCells,j-n2*numSubCells,k-n3*numSubCells);
                            val translation = imageTranslations(n1,n2,n3);
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
        
        timer.stop(TIMER_INDEX_DIRECT);
        return directEnergy();
    }

    /**
     * Generate an array of PointCharges, randomly distributed within a
     * size^3 cube centred at the origin.
     */
    private def generateAtoms(numAtoms:Long, seed:Int) : Rail[PointCharge] {
        val rand = seed > 0 ? new Random(seed) : new Random();
        Console.OUT.println("size of cluster =  " + size);
        val atoms = new Rail[PointCharge](numAtoms, (i:Long) => new PointCharge(Point3d((rand.nextDouble() - 0.5) * size, 
                                                    (rand.nextDouble() - 0.5) * size, 
                                                    (rand.nextDouble() - 0.5) * size), 
                                                i%2==0 ? 1.0 : -1.0));
        return atoms;
    }

    public static def main(args:Rail[String]) {
        val size:Double = 80.0;
        var numAtoms:Long;
        var cutoff:Double = 10.0;
        var randomSeed:Int = 0n;
        if (args.size > 0) {
            numAtoms = Long.parseLong(args(0));
            if (args.size > 1) {
                cutoff = Double.parseDouble(args(1));
                if (args.size > 2) {
                    randomSeed = Int.parseInt(args(2));
                }
            }
        } else {
            Console.ERR.println("usage: cutoff numAtoms [cutoff] [randomSeed]");
            return;
        }
        Console.OUT.println("electrostatics for " + numAtoms + " atoms with cutoff " + cutoff + " (total size " + size + ")");
        new ElectrostaticCutoff(size, cutoff).test(numAtoms, randomSeed);
    }

    public def test(numAtoms:Long, randomSeed:Int) {
        val atoms = generateAtoms(numAtoms, randomSeed);
        assignAtomsToSubCells(atoms);
        Console.OUT.println("Total PE: " + getEnergy());
        Console.OUT.printf("setup: %g seconds\n", (timer.mean(ElectrostaticCutoff.TIMER_INDEX_SETUP) as Double) / 1e9);
        Console.OUT.printf("getEnergy: %g seconds\n", (timer.mean(ElectrostaticCutoff.TIMER_INDEX_DIRECT) as Double) / 1e9);
    }

}
