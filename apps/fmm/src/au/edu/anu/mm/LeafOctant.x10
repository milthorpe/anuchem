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

import x10.util.ArrayList;
import x10.util.Pair;

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

    public def getAtomCharges():Rail[PointCharge] {
        val charges = new Array[PointCharge](atoms.size());
        for (i in 0..(atoms.size()-1)) {
            val atom = atoms(i);
            charges(i) = PointCharge(atom.centre, atom.charge);
        }
        return charges;
    }

    public def compareTo(b:LeafOctant):Int = id.compareTo(b.id);

    public def countOctants():Int = 1;

    public def ghostOctants():Int = 0;

    /**
     * Get the multipole representation of this octant, which is
     * simply the sum of the contributions of the particles in the octant.
     * N.B. must only be called once per pass
     */
    protected def upward(localData:PlaceLocalHandle[FmmLocalData], size:Double, dMax:UByte):Pair[Int,MultipoleExpansion] {
        //Console.OUT.println("at " + here + " LeafOctant.upward for " + id + " numAtoms = " + numAtoms);
        if (atoms.size() > 0) {
            multipoleExp.clear();
            localExp.clear();

            val p = multipoleExp.p;
            val centre = getCentre(size);
            for (i in 0..(atoms.size()-1)) {
                val atom = atoms(i);
                val atomLocation = centre.vector(atom.centre);
                // only one thread per octant, so unsafe addOlm is OK
                multipoleExp.addOlm(atom.charge, atomLocation, p);
            }

            atomic this.multipoleReady = true;

            sendMultipole(localData, dMax);

            return Pair[Int,MultipoleExpansion](atoms.size(), this.multipoleExp);
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
    }

    protected def downward(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, dMax:UByte):Double {
        //Console.OUT.println("at " + here + " LeafOctant.downward for " + id + " numAtoms = " + numAtoms);
        this.multipoleReady = false; // reset

        if (atoms.size() > 0) {
            constructLocalExpansion(localData, size, parentLocalExpansion);

            return getPotential(size, localData().locallyEssentialTree, dMax);
        } else {
            return 0.0;
        }        
    }

    /**
     * Returns the far-field potential of all charges within this octant due to
     * all charges in well-separated octant.
     * Updates forces on each particle due to long-range interactions.
     * @param size the side length of the full simulation box
     */
    private def getPotential(size:Double, myLET:LET, dMax:UByte):Double {
        val boxCentre = getCentre(size);

        var potential:Double = 0.0;
        for (atomIndex in 0..(atoms.size()-1)) {
            val atom = atoms(atomIndex);
            atom.force = Vector3d.NULL; // reset
            val locationWithinBox = atom.centre.vector(boxCentre);
            potential += localExp.calculatePotentialAndForces(atom, locationWithinBox);
        }

        potential += calculateDirectPotentialAndForces(size, myLET, dMax);

        return potential;
    }

    /** 
     * Calculates the potential and forces due to atoms in this octant and in
     * all non-well-separated octants.
     * @return the potential due to direct interactions
     */
    private def calculateDirectPotentialAndForces(size:Double, myLET:LET, dMax:UByte):Double {
        var directEnergy:Double = 0.0;
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

        val periodic = false; // TODO
        // direct calculation with all atoms in non-well-separated octants
        if (periodic) {
            val lowestLevelDim = Math.pow2(dMax);
            for (p in 0..(uList.size-1)) {
                val octantIndex2 = uList(p);
                val octant2Atoms = myLET.getAtomsForOctant(uList(p).getMortonId());
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
        } else {
            for (p in 0..(uList.size-1)) {
                val octant2Atoms = myLET.getAtomsForOctant(uList(p).getMortonId());
                if (octant2Atoms != null) {
                    //Console.OUT.println("at " + here + " calculating direct for " + id + " against " + uList(p));
                    for (octant2AtomsIndex in 0..(octant2Atoms.size-1)) {
                        val atom2 = octant2Atoms(octant2AtomsIndex);
                        for (atomIndex1 in 0..(atoms.size()-1)) {
                            val atom1 = atoms(atomIndex1);
                            val rVec = atom2.centre - atom1.centre;
                            val invR2 = 1.0 / rVec.lengthSquared();
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
     * Creates the U-list for this octant.
     * The U-list consists of all leaf octants not well-separated from this octant.
     */
    public def createUList(ws:Int) {
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
                        uList(i++) = OctantId(x2,y2,z2,id.level);
                    }
                }
            }
        }
    }

    public def getUList() = this.uList;

    public def toString(): String {
        return "LeafOctant " + id;
    }
}

