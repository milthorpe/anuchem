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

/**
 * This class represents a leaf node (with no children)
 * in the 3D division of space for the fast multipole method.
 * @author milthorpe
 */
public class FmmLeafBox extends FmmBox {
    private var atoms : Rail[PointCharge];

    /** The U-list consists of all leaf boxes not well-separated to this box. */
    private var uList : Rail[Point(3)];

    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : GlobalRef[FmmBox]) { 
        super(level, x, y, z, numTerms, parent);
    }

    public def getAtoms() = atoms;

    public def setAtoms(atoms : Rail[PointCharge]) {
        this.atoms = atoms;
    }

    protected def downward(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)]):Double {
        constructLocalExpansion(size, fmmOperators, parentLocalExpansion, locallyEssentialTree);

        return getPotential(size);
        
    }

    /**
     * Returns the far-field potential of all charges within this box due to
     * all charges in well-separated boxes.
     * @param size the side length of the full simulation box
     */
    private def getPotential(size:Double) : Double {
        val boxAtoms = getAtoms();
        val boxCentre = getCentre(size);
        val p = localExp.p;

        val chargeExpansion = new MultipoleExpansion(p);
        for (atomIndex in 0..(boxAtoms.size-1)) {
            val atom = boxAtoms(atomIndex);
            val locationWithinBox = atom.centre.vector(boxCentre);
            chargeExpansion.addOlm(atom.charge, locationWithinBox, p);
        }

        var potential : Double = 0.0;
        // TODO use lift/reduction?
        // TODO should be just:  for ([j,k] in terms.region) {
        for (j in 0..p) {
            for (k in -j..j) {
                potential += (localExp.terms(j,k) * chargeExpansion.terms(j,k)).re;
            }
        }
        return potential;
//        return terms.mapReduce[Complex,Double,Double](chargeExpansion.terms, 
//                                (a:Complex,b:Complex)=>(a*b).re, 
//                                (a:Double, b:Double)=>a+b, 0.0);
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
        // interact with "left half" of uList i.e. only boxes with x<=box.x
        val uList = new ArrayList[Point(3)]();
        for (x in Math.max(0,this.x-ws)..this.x) {
            for (y in Math.max(0,this.y-ws)..Math.min(levelDim-1,this.y+ws)) {
                for (z in Math.max(0,this.z-ws)..Math.min(levelDim-1,this.z+ws)) {
                    if (x < this.x || (x == this.x && y < this.y) || (x == this.x && y == this.y && z < this.z)) {
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
        // interact with "left half" of uList i.e. only boxes with x<=box.x
        val uList = new ArrayList[Point(3)]();
        for (x in (this.x-ws)..this.x) {
            for (y in (this.y-ws)..(this.y+ws)) {
                for (z in (this.z-ws)..(this.z+ws)) {
                    if (x < this.x || (x == this.x && y < this.y) || (x == this.x && y == this.y && z < this.z)) {
                        uList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        this.uList = uList.toArray();
    }
}

