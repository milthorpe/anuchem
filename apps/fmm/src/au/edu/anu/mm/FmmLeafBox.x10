/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.mm;

import x10.util.ArrayList;

import x10x.vector.Point3d;
import au.edu.anu.chem.PointCharge;
import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    private var atoms : Rail[MMAtom];

    /** The U-list consists of all leaf boxes not well-separated from this box. */
    private var uList : Rail[Point(3)];

    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : GlobalRef[FmmBox]) { 
        super(level, x, y, z, numTerms, parent);
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

    public def setAtoms(atoms : Rail[MMAtom]) {
        this.atoms = atoms;
    }

    /**
     * Get the multipole representation of this box, which as a leaf box is
     * simply the sum of the contributions of the particles in the box.
     * N.B. must only be called once per pass
     */
    protected def upward(size:Double, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)], periodic:Boolean) {
        val p = multipoleExp.p;
        val boxCentre = getCentre(size);
        val boxAtoms = getAtoms();
        multipoleExp.terms.clear();
        for (i in 0..(boxAtoms.size-1)) {
            val atom = boxAtoms(i);
            val atomLocation = boxCentre.vector(atom.centre);
            // only one thread per box, so unsafe addOlm is OK
            multipoleExp.addOlm(atom.charge, atomLocation, p);
        }

        sendMultipole(locallyEssentialTree, boxes, periodic);

        return multipoleExp;
    }

    protected def downward(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)]):Double {
        constructLocalExpansion(size, fmmOperators, parentLocalExpansion, locallyEssentialTree);
        val myLET = locallyEssentialTree();

        return getPotential(size, myLET);
        
    }

    /**
     * Returns the far-field potential of all charges within this box due to
     * all charges in well-separated boxes.
     * Updates forces on each particle due to long-range interactions.
     * @param size the side length of the full simulation box
     */
    private def getPotential(size:Double, myLET:LocallyEssentialTree) : Double {
        val boxAtoms = getAtoms();
        val boxCentre = getCentre(size);
        val p = localExp.p;

        val vExp = new MultipoleExpansion(p);
        val chargeExpansion = new MultipoleExpansion(p);
        for (atomIndex in 0..(boxAtoms.size-1)) {
            val atom = boxAtoms(atomIndex);
            val locationWithinBox = atom.centre.vector(boxCentre);
            val force = vExp.addOlmWithGradient(atom.charge, locationWithinBox, p, localExp);
            atom.force += force;
        }
        
        var potential:Double = 0.0;
        // TODO use lift/reduction?
        // TODO should be just:  for ([j,k] in terms.region) {
        for (j in 0..p) {
            potential += (localExp.terms(j,0) * vExp.terms(j,0)).re;
            for (k in 1..j) {
                val l_jk = localExp.terms(j,k);
                val v_jk = vExp.terms(j,k);
                potential += 2.0*(l_jk.re * v_jk.re) - 2.0*(l_jk.im * v_jk.im); // avoid conjugate/sign for mirror terms m < 0
            }
        }

        potential += calculateDirectPotentialAndForces(myLET);

        return potential;
//        return terms.mapReduce[Complex,Double,Double](chargeExpansion.terms, 
//                                (a:Complex,b:Complex)=>(a*b).re, 
//                                (a:Double, b:Double)=>a+b, 0.0);
    }

    /** 
     * Calculates the potential and forces due to atoms in this box and in
     * all non-well-separated boxes.
     * @return the potential due to direct interactions
     */
    private def calculateDirectPotentialAndForces(myLET:LocallyEssentialTree):Double {
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
        for (p in 0..(uList.size-1)) {
            val boxIndex2 = uList(p);
            // TODO - should be able to detect Point rank and inline
            val x2 = boxIndex2(0);
            val y2 = boxIndex2(1);
            val z2 = boxIndex2(2);
            val box2Atoms = cachedAtoms(x2, y2, z2);
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
        return directEnergy;
    }

    public def getUList() = this.uList;

    public def setUList(uList : Rail[Point(3)]) {
        this.uList = uList;
    }

    /**
     * Creates the U-list for this box.
     * The U-list consists of all leaf boxes not well-separated from this box.
     */
    public def createUList(ws : Int) {
        val levelDim = Math.pow2(this.level);
        val uList = new ArrayList[Point(3)]();
        for (x in Math.max(0,this.x-ws)..Math.min(levelDim-1,this.x+ws)) {
            for (y in Math.max(0,this.y-ws)..Math.min(levelDim-1,this.y+ws)) {
                for (z in Math.max(0,this.z-ws)..Math.min(levelDim-1,this.z+ws)) {
                    if (!(x==this.x && y==this.y && z==this.z)) {
                        uList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        this.uList = uList.toArray();
    }

    /**
     * Creates the U-list for this box for
     * use with the periodic FMM.
     * The U-list consists of all leaf boxes not well-separated from this box.
     */
    public def createUListPeriodic(ws : Int) {
        val levelDim = Math.pow2(this.level);
        val uList = new ArrayList[Point(3)]();
        for (x in (this.x-ws)..(this.x+ws)) {
            for (y in (this.y-ws)..(this.y+ws)) {
                for (z in (this.z-ws)..(this.z+ws)) {
                    if (!(x==this.x && y==this.y && z==this.z)) {
                        uList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        this.uList = uList.toArray();
    }
}

