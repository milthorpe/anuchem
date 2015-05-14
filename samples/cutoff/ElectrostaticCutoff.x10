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

import x10.array.Array_3;
import x10.regionarray.Array;
import x10.regionarray.Region;
import x10.util.ArrayList;
import x10.util.Random;
import x10.util.Timer;

/**
 * This class implements electrostatic potential calculation with cutoff
 * and periodic boundary conditions.
 */
public class ElectrostaticCutoff {
    public static class Atom(charge:Double) {
        public def this(charge:Double, x:Double, y:Double, z:Double) {
            property(charge);
            this.x = x;
            this.y = y;
            this.z = z;
        }

        // position
        public var x:Double;
        public var y:Double;
        public var z:Double;
        // velocity
        public var dx:Double;
        public var dy:Double;
        public var dz:Double;
        // force
        public var fx:Double;
        public var fy:Double;
        public var fz:Double;
    }

    /** The side length of the cubic simulation space. */
    private val size:Double;

    /** The direct sum cutoff distance */
    private val cutoff:Double;

    /** 
     * An array of box divisions within the unit cell, with a side length
     * equal to half the direct sum cutoff distance.  (N.B. if the unit cell
     * side length is not an exact multiple of the subcell side length, the
     * last box in each dimension will be smaller than the cutoff distance, 
     * resulting in anisotropy in the direct potential.)
     * Direct sums are only calculated between atoms in the same box and
     * the 26 neighbouring boxes.
     * Dimensions of the array region are (x,y,z)
     * TODO assumes cubic unit cell
     */
    private val subCells:Array_3[Rail[Atom]];

    /** The number of sub-cells per side of the unit cell. */
    private val numSubCells:Long;

    /**
     * @param size the length of a side of the cubic simulation space
     * @param cutoff the distance in Angstroms beyond which direct interactions are ignored
     */
    public def this(size:Double, cutoff : Double) {
        this.size = size;
        this.cutoff = cutoff;

        if (cutoff > size) {
            throw new IllegalArgumentException("error: cutoff " + cutoff + " is greater than size " + size);
        }
        if (size % (cutoff/2.0) != 0.0) {
            Console.ERR.println("warning: edge length " + size + " is not an exact multiple of (cutoff/2.0) " + (cutoff/2.0));
        }
        val numSubCells = Math.ceil(size / (cutoff/2.0)) as Int;
        val subCells = new Array_3[Rail[Atom]](numSubCells, numSubCells, numSubCells);
        this.subCells = subCells;
        this.numSubCells = numSubCells;
    }

    /**
     * Divide the atoms into a grid of sub-cells for direct sum calculation.
     * Each sub-cell is half the cutoff distance on every side.
     */
    private def assignAtomsToSubCells(atoms:Rail[Atom]) {
        val halfCutoff = (cutoff / 2.0);
        val subCellsTemp = new Array_3[ArrayList[Atom]](subCells.numElems_1, subCells.numElems_2, subCells.numElems_3, (i:Long,j:Long,k:Long) => new ArrayList[Atom]());
        val halfNumSubCells = this.numSubCells / 2;
        for (l in 0..(atoms.size-1)) {
            val atom = atoms(l);
            // get subcell i,j,k
            val i = (atom.x / halfCutoff) as Long + halfNumSubCells;
            val j = (atom.y / halfCutoff) as Long + halfNumSubCells;
            val k = (atom.z / halfCutoff) as Long + halfNumSubCells;
            subCellsTemp(i,j,k).add(atom);
        }
        for([i,j,k] in subCells.indices()) {
            subCells(i,j,k) = subCellsTemp(i,j,k).toRail();
        }
    }

    public def getEnergy() : Double {
        val cutoffSquared = cutoff*cutoff;
        val directEnergy = new Accumulator[Double](Reducible.SumReducer[Double]());
        finish for ([x,y,z] in subCells.indices()) async {
            val thisCell = subCells(x,y,z);
            var myDirectEnergy:Double = 0.0;
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
                            for (otherAtomIndex in 0..(otherCell.size-1)) {
                                val otherAtom = otherCell(otherAtomIndex);
                                val xj = otherAtom.x + n1*size;
                                val yj = otherAtom.y + n2*size;
                                val zj = otherAtom.z + n3*size;
                                val qj = otherAtom.charge;

                                val otherAtomCharge = otherAtom.charge;
                                for (thisAtomIndex in 0..(thisCell.size-1)) {
                                    val thisAtom = thisCell(thisAtomIndex);

                                    val rx = xj - thisAtom.x;
                                    val ry = yj - thisAtom.y;
                                    val rz = zj - thisAtom.z;
                                    
                                    val r2 = (rx*rx + ry*ry + rz*rz);
                                    if (r2 < cutoffSquared) {
                                        val invR:Double;
                                        val invR2:Double;
                                        invR2 = 1.0 / r2;
                                        invR = Math.sqrt(invR2);
                                        val qq = thisAtom.charge * qj;
                                        val e = invR * qq;
                                        myDirectEnergy += e;

                                        val forceScaling = e * invR2;
                                        thisAtom.fx -= forceScaling * rx;
                                        thisAtom.fy -= forceScaling * ry;
                                        thisAtom.fz -= forceScaling * rz;
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
                val qi = thisAtom.charge;
                val xi = thisAtom.x;
                val yi = thisAtom.y;
                val zi = thisAtom.z;
                var fix:Double = 0.0;
                var fiy:Double = 0.0;
                var fiz:Double = 0.0;

                for (j in 0..(i-1)) {
                    val otherAtom = thisCell(j);
                    val rx = otherAtom.x - xi;
                    val ry = otherAtom.y - yi;
                    val rz = otherAtom.z - zi;

                    val r2 = (rx*rx + ry*ry + rz*rz);

                    if (r2 < cutoffSquared) {
                        val invR:Double;
                        val invR2:Double;
                        invR2 = 1.0 / r2;
                        invR = Math.sqrt(invR2);
                        val qq = qi * otherAtom.charge;
                        val e = invR * qq;
                        myDirectEnergy += e;

                        val forceScaling = e * invR2;
                        fix -= forceScaling * rx;
                        fiy -= forceScaling * ry;
                        fiz -= forceScaling * rz;
                    }
                }
            }
            directEnergy <- myDirectEnergy;
        }
        
        return directEnergy();
    }

    /**
     * Generate atoms randomly distributed within a
     * size^3 cube centred at the origin.
     */
    private def generateAtoms(numAtoms:Long):Rail[Atom] {
        val rand = new Random();
        Console.OUT.println("size of cluster =  " + size);
        val randCoord = () => (rand.nextDouble() - 0.5) * size;
        val atoms = new Rail[Atom](numAtoms, 
            (i:Long) => new Atom(i%2==0 ? 1.0 : -1.0, randCoord(), randCoord(), randCoord()) );
        return atoms;
    }

    public static def main(args:Rail[String]) {
        val size:Double = 80.0;
        var numAtoms:Long;
        var cutoff:Double = 10.0;
        if (args.size > 0) {
            numAtoms = Long.parseLong(args(0));
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

    public def test(numAtoms:Long) {
        val setupStart = Timer.milliTime();
        val atoms = generateAtoms(numAtoms);
        assignAtomsToSubCells(atoms);
        val setupTime = Timer.milliTime() - setupStart;
        val energyStart = Timer.milliTime();
        val energy = getEnergy();
        val energyTime = Timer.milliTime() - energyStart;
        
        Console.OUT.println("Total PE: " + energy);
        Console.OUT.printf("setup: %d ms getEnergy: %d ms\n", setupTime, energyTime);
    }

}
