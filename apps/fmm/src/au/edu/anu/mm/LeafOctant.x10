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

import x10.compiler.Ifdef;
import x10.compiler.Ifndef;
import x10.compiler.Inline;
import x10.compiler.Native;
import x10.compiler.NativeCPPInclude;
import x10.util.ArrayList;
import x10.util.Pair;
import x10.util.Random;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.PointCharge;
import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
@NativeCPPInclude("bg_math.h")
public class LeafOctant extends Octant implements Comparable[LeafOctant] {
    public var atoms:ArrayList[MMAtom];
    private var sources:Rail[Double];

    /** The U-list consists of all leaf octants not well-separated from this octant. */
    private var uList:UList;

    public def this(id:OctantId, numTerms:Int, ws:Int, dMax:UByte) {
        super(id, numTerms, ws, dMax);
    }

    public def makeSources() {
        this.sources = new Rail[Double](atoms.size()*4);
        for (i in 0..(atoms.size()-1)) {
            val atom = atoms(i);
            val idx = i*4;
            sources(idx)   = atom.centre.i;
            sources(idx+1) = atom.centre.j;
            sources(idx+2) = atom.centre.k;
            sources(idx+3) = atom.charge;
        }
    }

    public def numAtoms() = atoms.size();

    public def getSources():Rail[Double] {
        return sources;
    }

    public def compareTo(b:LeafOctant):Int = id.compareTo(b.id);

    public def countOctants():Int = 1;

    public def ghostOctants():Int = 0;

    /**
     * Get the multipole representation of this octant, which is
     * simply the sum of the contributions of the particles in the octant.
     * N.B. must only be called once per pass
     */
    protected def upward():Pair[Long,MultipoleExpansion] {
        if (atoms.size() > 0) {
            multipoleExp.clear();
            localExp.clear();

            val size = FastMultipoleMethod.localData.size;
            val centre = getCentre(size);
            val plm = FmmScratch.getWorkerLocal().plm;
            for (i in 0..(atoms.size()-1)) {
                val atom = atoms(i);
                val atomLocation = centre.vector(atom.centre);
                // only one thread per octant, so unsafe addOlm is OK
                multipoleExp.addOlm(atom.charge, atomLocation, plm);
            }

            atomic this.multipoleReady = true;

            sendMultipole();

            return Pair[Long,MultipoleExpansion](atoms.size(), this.multipoleExp);
        } else {
            return Pair[Long,MultipoleExpansion](0L, null);
        }
    }

    protected def downward(parentLocalExpansion:LocalExpansion):Double {
        this.multipoleReady = false; // reset

        if (atoms.size() > 0) {
            addParentExpansion(parentLocalExpansion);

            val local = FastMultipoleMethod.localData;
            var potential: Double = farField(local.size);
@Ifndef("__EXCLUDE_NEAR__") {
            potential += nearField(local.size, local.locallyEssentialTree, local.dMax);
}

            return potential;
        } else {
            return 0.0;
        }  
    }

    /**
     * Calculates the potential and forces on atoms in this octant, due to 
     * all charges in well-separated octants, using the previously calculated
     * local expansion for this octant.
     * @param size the side length of the full simulation box
     * @return the potential due to far-field interactions
     */
    public def farField(size:Double):Double {
        val boxCentre = getCentre(size);

        val plm = FmmScratch.getWorkerLocal().plm;
        var potential:Double = 0.0;
        for (atomIndex in 0..(atoms.size()-1)) {
            val atom = atoms(atomIndex);
            atom.force = Vector3d.NULL; // reset
            val locationWithinBox = atom.centre.vector(boxCentre);
            potential += localExp.calculatePotentialAndForces(atom, locationWithinBox, plm);
        }

        return potential;
    }

    @Native("c++", "rsqrt(#a)")
    private @Inline native static def rsqrt(a:Double):Double;

    /** 
     * Calculates the potential and forces on atoms in this octant, due to 
     * all other atoms in this octant and in non-well-separated octants.
     * @param size the side length of the full simulation box
     * @param myLET the locally essential tree at this place
     * @param dMax the maximum number of levels in the tree
     * @return the potential due to near-field interactions
     */
    public def nearField(size:Double, myLET:LET, dMax:UByte):Double {
        var directEnergy:Double = 0.0;

        val periodic = false; // TODO
        // direct calculation with all atoms in non-well-separated octants
        for (mortonId in uList) {
            val oct2Data = myLET.getAtomDataForOctant(mortonId);
            if (oct2Data != null) {
                directEnergy += p2pKernel(oct2Data);
            }
        }

        for (atomIndex1 in 0..(atoms.size()-1)) {
            // direct calculation between all atoms in this octant
            val atom1 = atoms(atomIndex1);
            for (sameBoxAtomIndex in 0..(atomIndex1-1)) {
                val sameBoxAtom = atoms(sameBoxAtomIndex);
                val rVec = sameBoxAtom.centre - atom1.centre;
                val invR:Double;
                val invR2:Double;
@Ifdef("__BG__") {
                                invR = rsqrt(rVec.lengthSquared());
                                invR2 = invR * invR;
}
@Ifndef("__BG__") {
                                invR2 = 1.0 / rVec.lengthSquared();
                                invR = Math.sqrt(invR2);
}
                val e = atom1.charge * sameBoxAtom.charge * invR;
                directEnergy += 2.0 * e;
                val pairForce = e * invR2 * rVec;
                atom1.force += pairForce;
                sameBoxAtom.force -= pairForce;
            }
        }

        return directEnergy;
    }


    /**
     * Calculates forces and potential on atoms in this box due to atoms
     * in a neighbouring octant, where <code>oct2Data</code> is laid out 
     * sequentially in memory as a Rail[Double](N*4) repeating:
     * - atom x coord
     * - atom y coord
     * - atom z coord
     * - atom charge
     */
    private @Inline def p2pKernel(oct2Data:Rail[Double]) {
        var directEnergy:Double=0.0;
        for (i in 0..(atoms.size()-1)) {
            val atomI = atoms(i);
            val ci = atomI.centre;
            val xi = ci.i;
            val yi = ci.j;
            val zi = ci.k;
            val qi = atomI.charge;
            var fix:Double = atomI.force.i;
            var fiy:Double = atomI.force.j;
            var fiz:Double = atomI.force.k;

            for (var j:Long=0; j<oct2Data.size; j+=4) {
                val xj = oct2Data(j);
                val yj = oct2Data(j+1);
                val zj = oct2Data(j+2);
                val qj = oct2Data(j+3);

                val dx = xj-xi;
                val dy = yj-yi;
                val dz = zj-zi;
                val r2 = (dx*dx + dy*dy + dz*dz);
                val invR:Double;
                val invR2:Double;
@Ifdef("__BG__") {
                                invR = rsqrt(r2);
                                invR2 = invR * invR;
}
@Ifndef("__BG__") {
                                invR2 = 1.0 / r2;
                                invR = Math.sqrt(invR2);
}
                val qq = qi * qj;
                val e = invR * qq;
                directEnergy += e;

                val forceScaling = e * invR2;
                fix += forceScaling * dx;
                fiy += forceScaling * dy;
                fiz += forceScaling * dz;
            }
            atomI.force = Vector3d(fix, fiy, fiz);
        }
        return directEnergy;
    }

    /**
     * Gets the atom centre translation vector due to a lowest-level octant 
     * being in a neighbouring image, rather than the central unit cell.
     */
    private static def getTranslation(lowestLevelDim:Int, size:Double, x:Int, y:Int, z:Int) : Vector3d {
        var translationX:Double = 0.0;
        if (x >= lowestLevelDim) {
            translationX = size;
        } else if (x < 0) {
            translationX = -size;
        }

        var translationY:Double = 0;
        if (y >= lowestLevelDim) {
            translationY = size;
        } else if (y < 0) {
            translationY = -size;
        }

        var translationZ:Double = 0;
        if (z >= lowestLevelDim) {
            translationZ = size;
        } else if (z < 0) {
            translationZ = -size;
        }
        return Vector3d(translationX, translationY, translationZ);
    }

    /** 
     * Estimates the number of octants in a leaf octant's U-list,
     * given the octant id. 
     */
    public static def estimateUListSize(mortonId:UInt, dMax:UByte):Int {
        val id = OctantId.getFromMortonId(mortonId);
        return estimateUListSize(id, dMax);
    }

    public static def estimateUListSize(id:OctantId, dMax:UByte):Int {
        val maxExtent = (1U << dMax) - 1U;
        val xExtent = (id.x > 0U && id.x < maxExtent) ? 3 : 2;
        val yExtent = (id.y > 0U && id.y < maxExtent) ? 3 : 2;
        val zExtent = (id.z > 0U && id.z < maxExtent) ? 3 : 2;
        return xExtent * yExtent * zExtent;
    }

    /**
     * Returns a cost estimate per interaction (in ns) of U-list calculation.
     * @param q the average density per lowest level box
     */
    public static def estimateUListCost(q:Int):Long {
        val dummyOctant = new LeafOctant(new OctantId(0UY, 0UY, 0UY, 0UY), 1, 1, 0UY);
        dummyOctant.atoms = new ArrayList[MMAtom](q);
        val rand = new Random();
        for (i in 0..(q-1)) {
            // create dummy singly-charged ions in random positions offset by [40,40,40]
            val atom = new MMAtom(new Point3d(rand.nextDouble()+40.0, rand.nextDouble()+40.0, rand.nextDouble()+40.0), 1.0, 1.0);
            dummyOctant.atoms(i) = atom;
        }

        // interact with dummy singly-charged ions in random positions offset by [39,39,39]
        val dummyData = new Rail[Double](q*4, 
                (i:Long) => (i%4==3L)? 1.0 : rand.nextDouble()+39.0);

        val start = System.nanoTime();
        for (i in 0..1000) {
            dummyOctant.p2pKernel(dummyData);
        }
        val stop = System.nanoTime();

        val interactions = 1000 * q * q;
        val perInt = (stop-start) / interactions;
        //Console.OUT.println("for " + id + " q = " + atoms.size() + " interactions = " + interactions + " perInt = " + perInt);
        return perInt;
    }

    public def createUList(ws:Int) {
        uList = new UList(id, ws);
    }

    public def getUList() = uList;

    /** 
     * The U-list consists of all leaf octants not well-separated from this octant.
     */
    public class UList implements Iterable[UInt] {
        val level:UByte;
        val minX:UByte;
        val maxX:UByte;
        val minY:UByte;
        val maxY:UByte;
        val minZ:UByte;
        val maxZ:UByte;

        public def this(id:OctantId, ws:Int) {
            val levelDim = 1UY << id.level;
            minX = Math.max(0,id.x-ws) as UByte;
            maxX = Math.min(levelDim-1,id.x+ws) as UByte;
            minY = Math.max(0,id.y-ws) as UByte;
            maxY = Math.min(levelDim-1,id.y+ws) as UByte;
            minZ = Math.max(0,id.z-ws) as UByte;
            maxZ = Math.min(levelDim-1,id.z+ws) as UByte;
            this.level = id.level;
        }

        public def iterator() = new UListIterator();

        public class UListIterator implements Iterator[UInt] {
            var x:UByte;
            var y:UByte;
            var z:UByte;
            public def this() {
                x = minX;
                y = minY;
                z = minZ;
            }

            public def hasNext():Boolean {
                if (x <= maxX) {
                    if (!(x==id.x && y==id.y && z==id.z)) {
                        return true;
                    } else {
                        moveToNext();
                        if (x <= maxX && !(x==id.x && y==id.y && z==id.z)) return true;
                    }
                }
                return false;
            }
            
            public def next():UInt {
                if (x <= maxX && !(x==id.x && y==id.y && z==id.z)) {
                    val res = OctantId.getMortonId(x, y, z, level);
                    moveToNext();
                    return res;
                } else {
                    throw new UnsupportedOperationException("reached end of uList for " + LeafOctant.this.id);
                }
            }

            private def moveToNext() {
                do {
                    if (z < maxZ) {
                        z++;
                    } else if (y < maxY) {
                        z = minZ;
                        y++;
                    } else {
                        z = minZ;
                        y = minY;
                        x++;
                    }
                } while(x <= maxX && (x==id.x && y==id.y && z==id.z));
            }
        }
    }

    public def toString(): String {
        return "LeafOctant " + id;
    }
}

