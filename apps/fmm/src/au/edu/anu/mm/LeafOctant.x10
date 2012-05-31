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
    private var atoms:Rail[MMAtom];
    public var atomList:ArrayList[MMAtom];

    /** The U-list consists of all leaf boxes not well-separated from this box. */
    private var uList:Rail[OctantId];

    public def this(id:OctantId, numTerms:Int) {
        super(id, numTerms);
        atomList = new ArrayList[MMAtom]();
    }

    public def getAtomCharges():Rail[PointCharge] {
        val charges = new Array[PointCharge](atoms.size);
        for (i in atoms) {
            val atom = atoms(i);
            charges(i) = PointCharge(atom.centre, atom.charge);
        }
        return charges;
    }

    public def getAtoms() = atoms;

    public def setAtoms(atoms:Rail[MMAtom]) {
        this.atoms = atoms;
        numAtoms = atoms.size;
    }

    public def compareTo(b:LeafOctant):Int = id.compareTo(b.id);

    /**
     * Get the multipole representation of this octant, which is
     * simply the sum of the contributions of the particles in the octant.
     * N.B. must only be called once per pass
     */
    protected def upward(size:Double, fmmOperators:FmmOperators, locallyEssentialTree:LET, periodic:Boolean):Pair[Int,MultipoleExpansion] {
        if (atoms.size > 0) {
            val p = multipoleExp.p;
            val centre = getCentre(size);
            multipoleExp.terms.clear();
            for (i in 0..(atoms.size-1)) {
                val atom = atoms(i);
                val atomLocation = centre.vector(atom.centre);
                // only one thread per octant, so unsafe addOlm is OK
                multipoleExp.addOlm(atom.charge, atomLocation, p);
            }

            // TODO sendMultipole(locallyEssentialTree, boxes, periodic);

            return Pair[Int,MultipoleExpansion](numAtoms, this.multipoleExp);
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
    }

    protected def downward(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:FmmOperators, locallyEssentialTree:LET, numLevels:Int, periodic:Boolean):Double {
        if (numAtoms > 0) {
            constructLocalExpansion(size, fmmOperators, parentLocalExpansion, locallyEssentialTree);

            return getPotential(size, locallyEssentialTree, numLevels, periodic);
        } else {
            return 0.0;
        }        
    }

    /**
     * Returns the far-field potential of all charges within this box due to
     * all charges in well-separated boxes.
     * Updates forces on each particle due to long-range interactions.
     * @param size the side length of the full simulation box
     */
    private def getPotential(size:Double, myLET:LET, numLevels:Int, periodic:Boolean):Double {
        val boxCentre = getCentre(size);

        var potential:Double = 0.0;
        for (atomIndex in 0..(atoms.size-1)) {
            val atom = atoms(atomIndex);
            atom.force = Vector3d.NULL; // reset
            val locationWithinBox = atom.centre.vector(boxCentre);
            potential += localExp.calculatePotentialAndForces(atom, locationWithinBox);
        }

        potential += calculateDirectPotentialAndForces(size, myLET, numLevels, periodic);

        return potential;
    }

    /** 
     * Calculates the potential and forces due to atoms in this box and in
     * all non-well-separated boxes.
     * @return the potential due to direct interactions
     */
    private def calculateDirectPotentialAndForces(size:Double, myLET:LET, numLevels:Int, periodic:Boolean):Double {
        val cachedAtoms = myLET.cachedAtoms;
        var directEnergy:Double = 0.0;
        for (atomIndex1 in 0..(atoms.size-1)) {
            // direct calculation between all atoms in this box
            val atom1 = atoms(atomIndex1);
            for (sameBoxAtomIndex in 0..(atomIndex1-1)) {
                val sameBoxAtom = atoms(sameBoxAtomIndex);
                val rVec = sameBoxAtom.centre - atom1.centre;
                val r2 = rVec.lengthSquared();
                val r = Math.sqrt(r2);
                val pairEnergy = 2.0 * atom1.charge * sameBoxAtom.charge / r;
                directEnergy += pairEnergy;
                val pairForce = (atom1.charge * sameBoxAtom.charge / r2 / r) * rVec;
                atom1.force += pairForce;
                sameBoxAtom.force -= pairForce;
            }
        }

        // direct calculation with all atoms in non-well-separated boxes
        if (periodic) {
            val lowestLevelDim = Math.pow2(numLevels);
            for (p in 0..(uList.size-1)) {
                val box2Atoms = cachedAtoms(p);
                if (box2Atoms != null) {
                    val octantIndex2 = uList(p);
                    val translation = getTranslation(lowestLevelDim, size, octantIndex2.x, octantIndex2.y, octantIndex2.z);
                    for (otherBoxAtomIndex in 0..(box2Atoms.size-1)) {
                        val atom2 = box2Atoms(otherBoxAtomIndex);
                        val translatedCentre = atom2.centre + translation;
                        for (atomIndex1 in 0..(atoms.size-1)) {
                            val atom1 = atoms(atomIndex1);
                            val rVec = translatedCentre - atom1.centre;
                            val r2 = rVec.lengthSquared();
                            val r = Math.sqrt(r2);
                            if (r != 0.0) { // don't include dipole-balancing charges at same point
                                directEnergy += atom1.charge * atom2.charge / r;
                                val pairForce = (atom1.charge * atom2.charge / r2 / r) * rVec;
                                atom1.force += pairForce;
                            }
                        }
                    }
                }
            }
        } else {
            for (p in 0..(uList.size-1)) {
                val box2Atoms = cachedAtoms(p);
                if (box2Atoms != null) {
                    for (otherBoxAtomIndex in 0..(box2Atoms.size-1)) {
                        val atom2 = box2Atoms(otherBoxAtomIndex);
                        for (atomIndex1 in 0..(atoms.size-1)) {
                            val atom1 = atoms(atomIndex1);
                            val rVec = atom2.centre - atom1.centre;
                            val r2 = rVec.lengthSquared();
                            val r = Math.sqrt(r2);
                            directEnergy += atom1.charge * atom2.charge / r;
                            val pairForce = (atom1.charge * atom2.charge / r2 / r) * rVec;
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
     * Creates the U-list for this box.
     * The U-list consists of all leaf boxes not well-separated from this box.
     */
    public def createUList(ws:Int) {
        val levelDim = Math.pow2(id.level);
        val w = 2*ws+1;
        uList = new Array[OctantId](w*w*w-1);
        var i:Int = 0;
        for (x in Math.max(0,id.x-ws)..Math.min(levelDim-1,id.x+ws)) {
            val x2 = x as UShort;
            for (y in Math.max(0,id.y-ws)..Math.min(levelDim-1,id.y+ws)) {
                val y2 = y as UShort;
                for (z in Math.max(0,id.z-ws)..Math.min(levelDim-1,id.z+ws)) {
                    val z2 = z as UShort;
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

