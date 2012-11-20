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
    public static val TIMER_INDEX_SETUP : Int = 2;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(3);

    /** The direct sum cutoff distance in Angstroms */
    private val cutoff : Double;

    /** 
     * Translation vectors for neighbouring unit cells 
     * (the 26 cells surrounding the origin cell).
     * Dimensions are:
     * 0: x translation (difference between x-coordinate of sub-cells
     * 1: y translation
     * 2: z translation
     */
    private val imageTranslations : Array[Vector3d](3){rect};


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
    private val subCells : Array[Rail[PointCharge]](3){rect,zeroBased};

    /** The number of sub-cells per side of the unit cell. */
    private val numSubCells : Int;

    /**
     * @param size the length of a side of the cubic simulation space
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(size:Double,
            cutoff : Double) {

        this.cutoff = cutoff;
        this.imageTranslations = new Array[Vector3d](-1..1 * -1..1 * -1..1, 
                ([i,j,k] : Point(3)) => Vector3d(i*size, j*size, k*size));

        if (cutoff > size) {
            throw new IllegalArgumentException("error: cutoff " + cutoff + " is greater than size " + size);
        }
        if (size % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + size + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(size / (cutoff/2.0)) as Int;
        val subCellRegion = 0..(numSubCells-1) * 0..(numSubCells-1) * 0..(numSubCells-1);
        val subCells = new Array[Rail[PointCharge]](subCellRegion);
        this.subCells = subCells;
        this.numSubCells = numSubCells;
    }

    /**
     * Divide the atoms into a grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def assignAtomsToSubCells(atoms: Rail[MMAtom]) {
        timer.start(TIMER_INDEX_SETUP);
        val halfCutoff = (cutoff / 2.0);
        val subCellsTemp = new Array[ArrayList[PointCharge]](subCells.region, (Point) => new ArrayList[PointCharge]());
        val halfNumSubCells = this.numSubCells / 2;
        for (l in 0..(atoms.size-1)) {
            val atom = atoms(l);
            val centre = atom.centre;
            val charge = atom.charge;
            // get subcell i,j,k
            val i = (centre.i / halfCutoff) as Int + halfNumSubCells;
            val j = (centre.j / halfCutoff) as Int + halfNumSubCells;
            val k = (centre.k / halfCutoff) as Int + halfNumSubCells;
            subCellsTemp(i,j,k).add(new PointCharge(centre, charge));
        }
        for([i,j,k] in subCells) {
            subCells(i,j,k) = subCellsTemp(i,j,k).toArray();
        }
        timer.stop(TIMER_INDEX_SETUP);
    }

    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_DIRECT);
        val cutoffSquared = cutoff*cutoff;
        val directEnergy = new Accumulator[Double](Reducible.SumReducer[Double]());
        finish for ([x,y,z] in subCells) async {
            val thisCell = subCells(x,y,z);
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
        assignAtomsToSubCells(atoms(0));
        Console.OUT.println("Total PE: " + getEnergy());
        logTime("Setup",             ElectrostaticCutoff.TIMER_INDEX_SETUP,         timer);
        logTime("Direct",            ElectrostaticCutoff.TIMER_INDEX_DIRECT,        timer);
    }

}
