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
package au.edu.anu.mm;

import x10.compiler.Inline;
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
public class LeafOctant extends Octant implements Comparable[LeafOctant] {
    public var atoms:ArrayList[MMAtom];

    /** The U-list consists of all leaf octants not well-separated from this octant. */
    private var uList:Rail[OctantId];

    public def this(id:OctantId, numTerms:Int) {
        super(id, numTerms);
        atoms = new ArrayList[MMAtom]();
    }

    public def numAtoms() = atoms.size();

    public def getAtomData():Rail[Double] {
        val flat = new Array[Double](atoms.size()*4);
        for (i in 0..(atoms.size()-1)) {
            val atom = atoms(i);
            val idx = i*4;
            flat(idx)   = atom.centre.i;
            flat(idx+1) = atom.centre.j;
            flat(idx+2) = atom.centre.k;
            flat(idx+3) = atom.charge;
        }
        return flat;
    }

    public def compareTo(b:LeafOctant):Int = id.compareTo(b.id);

    public def countOctants():Int = 1;

    public def ghostOctants():Int = 0;

    /**
     * Get the multipole representation of this octant, which is
     * simply the sum of the contributions of the particles in the octant.
     * N.B. must only be called once per pass
     */
    protected def upward(localData:PlaceLocalHandle[FmmLocalData]):Pair[Int,MultipoleExpansion] {
        //Console.OUT.println("at " + here + " LeafOctant.upward for " + id + " numAtoms = " + numAtoms);
        if (atoms.size() > 0) {
            multipoleExp.clear();
            localExp.clear();

            val size = localData().size;
            val p = multipoleExp.p;
            val centre = getCentre(size);
            for (i in 0..(atoms.size()-1)) {
                val atom = atoms(i);
                val atomLocation = centre.vector(atom.centre);
                // only one thread per octant, so unsafe addOlm is OK
                multipoleExp.addOlm(atom.charge, atomLocation, p);
            }

            atomic this.multipoleReady = true;

            sendMultipole(localData);

            return Pair[Int,MultipoleExpansion](atoms.size(), this.multipoleExp);
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
    }

    protected def downward(localData:PlaceLocalHandle[FmmLocalData], parentLocalExpansion:LocalExpansion):Double {
        //Console.OUT.println("at " + here + " LeafOctant.downward for " + id + " numAtoms = " + numAtoms);
        this.multipoleReady = false; // reset

        if (atoms.size() > 0) {
            constructLocalExpansion(localData, parentLocalExpansion);

            val local = localData();
            var potential: Double = farField(local.size);
            potential += nearField(local.size, local.locallyEssentialTree, local.dMax);

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

        var potential:Double = 0.0;
        for (atomIndex in 0..(atoms.size()-1)) {
            val atom = atoms(atomIndex);
            atom.force = Vector3d.NULL; // reset
            val locationWithinBox = atom.centre.vector(boxCentre);
            potential += localExp.calculatePotentialAndForces(atom, locationWithinBox);
        }

        return potential;
    }

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
        if (periodic) {
/*
            val lowestLevelDim = Math.pow2(dMax);
            for (p in 0..(uList.size-1)) {
                val octantIndex2 = uList(p);
                val octant2Atoms = myLET.getAtomDataForOctant(uList(p).getMortonId());
                if (octant2Atoms != null) {
                    val translation = getTranslation(lowestLevelDim, size, octantIndex2.x, octantIndex2.y, octantIndex2.z);
                    for (octant2AtomsIndex in 0..(octant2Atoms.size-1)) {
                        val atom2 = octant2Atoms(octant2AtomsIndex);
                        val translatedCentre = atom2.centre + translation;
                        for (atomIndex1 in 0..(atoms.size()-1)) {
                            val atom1 = atoms(atomIndex1);
                            val rVec = translatedCentre - atom1.centre;
                            val r2 = rVec.lengthSquared();
                            if (r2 != 0.0) { // don't include dipole-balancing charges at same point
                                val invR2 = 1.0 / r2;
                                val invR = Math.sqrt(invR2);
                                val e = atom1.charge * atom2.charge * invR;
                                directEnergy += e;
                                val pairForce = e * invR2 * rVec;
                                atom1.force += pairForce;
                            }
                        }
                    }
                }
            }
*/
        } else {
            //Console.OUT.println(id + " uList " + uList.size);
            for (p in 0..(uList.size-1)) {
                val oct2Data = myLET.getAtomDataForOctant(uList(p).getMortonId());
                if (oct2Data != null) {
                    directEnergy += p2pKernel(oct2Data);
                }
            }
        }

        for (atomIndex1 in 0..(atoms.size()-1)) {
            // direct calculation between all atoms in this octant
            val atom1 = atoms(atomIndex1);
            for (sameBoxAtomIndex in 0..(atomIndex1-1)) {
                val sameBoxAtom = atoms(sameBoxAtomIndex);
                val rVec = sameBoxAtom.centre - atom1.centre;
                val invR2 = 1.0 / rVec.lengthSquared();
                val invR = Math.sqrt(invR2);
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

            //Console.OUT.println("at " + here + " calculating direct for " + id + " against " + uList(p));
            for (var j:Int=0; j<oct2Data.size; j+=4) {
                val xj = oct2Data(j);
                val yj = oct2Data(j+1);
                val zj = oct2Data(j+2);
                val qj = oct2Data(j+3);

                val dx = xj-xi;
                val dy = yj-yi;
                val dz = zj-zi;
                val r2 = (dx*dx + dy*dy + dz*dz);
                val invR2 = 1.0 / r2;
                val invR = Math.sqrt(invR2);
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
        val maxExtent = (1U << dMax) - 1U;
        val id = OctantId.getFromMortonId(mortonId);
        val xExtent = (id.x > 0U && id.x < maxExtent) ? 3 : 2;
        val yExtent = (id.y > 0U && id.y < maxExtent) ? 3 : 2;
        val zExtent = (id.z > 0U && id.z < maxExtent) ? 3 : 2;
        //Console.OUT.println("uList for " + id + " size = " + xExtent * yExtent * zExtent);
        return xExtent * yExtent * zExtent;
    }

    /**
     * Returns a cost estimate per interaction (in ns) of U-list calculation.
     */
    public def estimateUListCost():Long {
        // use number of particles in this box as an estimate for others
        val rand = new Random();
        // create dummy singly-charged ions in random positions offset by [40,40,40]
        val dummyData = new Rail[Double](atoms.size()*4, 
                (i:Int) => (i%4==3)? 1.0 : rand.nextDouble()+40.0);
        val start = System.nanoTime();
        for (i in 0..26) {
            p2pKernel(dummyData);
        }
        val stop = System.nanoTime();
        val interactions = 27 * atoms.size() * atoms.size();
        val perInt = (stop-start) / interactions;
        //Console.OUT.println("for " + id + " q = " + atoms.size() + " interactions = " + interactions + " perInt = " + perInt);
        return perInt;
    }

    /**
     * Creates the U-list for this octant.
     * The U-list consists of all leaf octants not well-separated from this octant.
     */
    public def createUList(ws:Int) {
        val tempUList = new ArrayList[OctantId]();
        val levelDim = Math.pow2(id.level);
        val w = 2*ws+1;
        uList = new Array[OctantId](w*w*w-1);
        var i:Int = 0;
        for (x in Math.max(0,id.x-ws)..Math.min(levelDim-1,id.x+ws)) {
            val x2 = x as UByte;
            for (y in Math.max(0,id.y-ws)..Math.min(levelDim-1,id.y+ws)) {
                val y2 = y as UByte;
                for (z in Math.max(0,id.z-ws)..Math.min(levelDim-1,id.z+ws)) {
                    val z2 = z as UByte;
                    if (!(x2==id.x && y2==id.y && z2==id.z)) {
                        tempUList.add(OctantId(x2,y2,z2,id.level));
                    }
                }
            }
        }
        this.uList = tempUList.toArray();
    }

    public def getUList() = this.uList;

    public def toString(): String {
        return "LeafOctant " + id;
    }
}

