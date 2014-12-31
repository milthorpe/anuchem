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
import au.edu.anu.chem.mm.AtomType;
import au.edu.anu.chem.mm.ParticleData;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
@NativeCPPInclude("bg_math.h")
public class LeafOctant extends Octant {
    public val atoms = new ArrayList[Long]();
    private var sources:Rail[Double];

    /** The U-list consists of all leaf octants not well-separated from this octant. */
    private var uList:UList;

    public def this(id:OctantId, numTerms:Int, ws:Int, dMax:UByte) {
        super(id, numTerms, ws, dMax);
    }

    public def makeSources(particleData:ParticleData) {
        this.sources = new Rail[Double](atoms.size()*4);
        for (i in 0..(atoms.size()-1)) {
            val atomIndex = atoms(i);
            val idx = i*4;
            sources(idx)   = particleData.x(atomIndex).i;
            sources(idx+1) = particleData.x(atomIndex).j;
            sources(idx+2) = particleData.x(atomIndex).k;
            sources(idx+3) = particleData.atomTypes(particleData.atomTypeIndex(atomIndex)).charge;
        }
    }

    public def numAtoms() = atoms.size();

    public def getSources():Rail[Double] {
        return sources;
    }

    public def countOctants():Int = 1n;

    public def ghostOctants():Int = 0n;

    /**
     * Get the multipole representation of this octant, which is
     * simply the sum of the contributions of the particles in the octant.
     * N.B. must only be called once per pass
     */
    protected def upward():Pair[Long,MultipoleExpansion] {
        if (atoms.size() > 0) {
            multipoleExp.clear();
            localExp.clear();

            val local = FastMultipoleMethod.localData;
            val size = local.size;
            val particleData = local.particleData;
            val centre = getCentre(size);
            val plm = FmmScratch.getWorkerLocal().plm;
            for (atomIdx in atoms) {
                val atomLocation = centre.vector(particleData.x(atomIdx));
                // only one thread per octant, so unsafe addOlm is OK
                val atomCharge = particleData.atomTypes(particleData.atomTypeIndex(atomIdx)).charge;
                multipoleExp.addOlm(atomCharge, atomLocation, plm);
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
            var potential: Double = farField(local.size, local.particleData);
@Ifndef("__EXCLUDE_NEAR__") {
            potential += nearField(local.size, local.particleData, local.locallyEssentialTree, local.dMax);
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
    public def farField(size:Double, particleData:ParticleData):Double {
        val boxCentre = getCentre(size);

        val plm = FmmScratch.getWorkerLocal().plm;
        var potential:Double = 0.0;
        for (atomIdx in atoms) {
            val atomCharge = particleData.atomTypes(particleData.atomTypeIndex(atomIdx)).charge;
            particleData.fx(atomIdx) = Vector3d.NULL; // reset
            val locationWithinBox = particleData.x(atomIdx).vector(boxCentre);
            potential += localExp.calculatePotentialAndForces(atomIdx, particleData, locationWithinBox, plm);
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
    public def nearField(size:Double, particleData:ParticleData, myLET:LET, dMax:UByte):Double {
        var directEnergy:Double = 0.0;

        val periodic = false; // TODO
        // direct calculation with all atoms in non-well-separated octants
        for (mortonId in uList) {
            val oct2Data = myLET.getAtomDataForOctant(mortonId);
            if (oct2Data != null) {
                directEnergy += p2pKernel(oct2Data, particleData);
            }
        }

        for (a1 in 0..(atoms.size()-1)) {
            val atomIndex1 = atoms(a1);
            // direct calculation between all atoms in this octant
            for (a2 in 0..(a1-1)) {
                val atomIndex2 = atoms(a2);
                val rVec = particleData.x(atomIndex2) - particleData.x(atomIndex1);
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
                val e = particleData.atomTypes(particleData.atomTypeIndex(atomIndex1)).charge
                      * particleData.atomTypes(particleData.atomTypeIndex(atomIndex2)).charge
                      * invR;
                directEnergy += 2.0 * e;
                val pairForce = e * invR2 * rVec;
                particleData.fx(atomIndex1) += pairForce;
                particleData.fx(atomIndex2) -= pairForce;
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
    private @Inline def p2pKernel(oct2Data:Rail[Double], particleData:ParticleData) {
        var directEnergy:Double=0.0;
        for (i in atoms) {
            val ci = particleData.x(i);
            val xi = ci.i;
            val yi = ci.j;
            val zi = ci.k;
            val qi = particleData.atomTypes(particleData.atomTypeIndex(i)).charge;
            val fi = particleData.fx(i);
            var fix:Double = fi.i;
            var fiy:Double = fi.j;
            var fiz:Double = fi.k;

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
            particleData.fx(i) = Vector3d(fix, fiy, fiz);
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
        val maxExtent = (1UN << dMax) - 1UN;
        val xExtent = (id.x > 0U && id.x < maxExtent) ? 3n : 2n;
        val yExtent = (id.y > 0U && id.y < maxExtent) ? 3n : 2n;
        val zExtent = (id.z > 0U && id.z < maxExtent) ? 3n : 2n;
        return xExtent * yExtent * zExtent;
    }

    /**
     * Returns a cost estimate per interaction (in ns) of U-list calculation.
     * @param q the average density per lowest level box
     */
    public static def estimateUListCost(q:Int):Long {
        val dummyOctant = new LeafOctant(new OctantId(0UY, 0UY, 0UY, 0UY), 1n, 1n, 0UY);
        val dummyParticles = new ParticleData();
        dummyParticles.atomTypes = [ AtomType("default", -1n, 1.0, 1.0) as AtomType ];
        val rand = new Random();
        for (i in 0..(q-1)) {
            // create dummy singly-charged ions in random positions offset by [40,40,40]
            dummyParticles.addAtom(i, 0n, new Point3d(rand.nextDouble()+40.0, rand.nextDouble()+40.0, rand.nextDouble()+40.0));
        }

        // interact with dummy singly-charged ions in random positions offset by [39,39,39]
        val dummyData = new Rail[Double](q*4, 
                (i:Long) => (i%4==3L)? 1.0 : rand.nextDouble()+39.0);

        val start = System.nanoTime();
        for (i in 0..1000) {
            dummyOctant.p2pKernel(dummyData, dummyParticles);
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
            minX = Math.max(0n,id.x-ws) as UByte;
            maxX = Math.min(levelDim-1n,id.x+ws) as UByte;
            minY = Math.max(0n,id.y-ws) as UByte;
            maxY = Math.min(levelDim-1n,id.y+ws) as UByte;
            minZ = Math.max(0n,id.z-ws) as UByte;
            maxZ = Math.min(levelDim-1n,id.z+ws) as UByte;
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

