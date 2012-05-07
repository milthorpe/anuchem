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
import x10.util.HashSet;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents a box in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public class FmmBox {
    public val parent : GlobalRef[FmmBox];

    public val level : Int;
    public val x : Int;
    public val y : Int;
    public val z : Int;

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private var vList : Rail[Point(3)];

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp : MultipoleExpansion;

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp : LocalExpansion;

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(level : Int, x : Int, y : Int, z : Int, numTerms : Int, parent : GlobalRef[FmmBox]) {
        this.level = level;
        this.x = x;
        this.y = y;
        this.z = z;
        this.parent = parent;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public def getCentre(size : Double) : Point3d {
        dim : Int = Math.pow2(level);
        sideLength : Double = size / dim;
        offset : Double = 0.5 * size;
        return Point3d( (x + 0.5) * sideLength - offset,
                        (y + 0.5) * sideLength - offset,
                        (z + 0.5) * sideLength - offset);
    }

    /**
     * Returns true if this box is well-separated from <code>x,y,z</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public def wellSeparated(ws : Int, x2 : Int, y2 : Int, z2 : Int) : Boolean {
        return Math.abs(x - x2) > ws 
            || Math.abs(y - y2) > ws 
            || Math.abs(z - z2) > ws;
    }

    /**
     * Returns true if this box is well-separated from <code>box2</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public def wellSeparated(ws : Int, box2 : FmmBox) : Boolean {
        return Math.abs(x - box2.x) > ws 
            || Math.abs(y - box2.y) > ws 
            || Math.abs(z - box2.z) > ws;
    }

    private static def upwardIfNotNull(size:Double, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)], childLevel:Int, x:Int, y:Int, z:Int, periodic:Boolean) : MultipoleExpansion {
        val childLevelBoxes = boxes(childLevel);
        val childBox = childLevelBoxes(x,y,z);
        if (childBox != null) {
            return childBox.upward(size, fmmOperators, locallyEssentialTree, boxes, periodic);
        } else {
            return null;
        }
    }

    /** 
     * For each non-leaf box, combines multipole expansions for <= 8 child
     * boxes into a single multipole expansion for the parent box.
     */
    protected def upward(size:Double, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)], periodic:Boolean) {
        val childLevel = level+1;
        val childLevelBoxes = boxes(childLevel);

        val childExpansions = new Array[MultipoleExpansion](8);
        finish {
            var i:Int=0;
            for (x2 in (2*x)..(2*x+1)) {
                for (y2 in (2*y)..(2*y+1)) {
                    for (z2 in (2*z)..(2*z+1)) {
                        val childHome = childLevelBoxes.dist(x2,y2,z2);
                        if (childHome == here) {
                            childExpansions(i++) = upwardIfNotNull(size, fmmOperators, locallyEssentialTree, boxes, childLevel, x2, y2, z2, periodic);
                        } else {
                            childExpansions(i++) = at (childHome) upwardIfNotNull(size, fmmOperators, locallyEssentialTree, boxes, childLevel, x2, y2, z2, periodic);
                        }
                    }
                }
            }
            multipoleExp.terms.clear();
        }

        val myOperators = fmmOperators();
        val myComplexK = myOperators.complexK;
        val myWignerA = myOperators.wignerA;
        val halfSideLength = size / Math.pow2(childLevel+1);
        val numTerms = multipoleExp.p;
        val scratch = new MultipoleExpansion(numTerms);    
        val scratch_array = new Array[Complex](numTerms+1);
        var i:Int=0;
        for (x2 in (2*x)..(2*x+1)) {
            for (y2 in (2*y)..(2*y+1)) {
                for (z2 in (2*z)..(2*z+1)) {
                    val childExp = childExpansions(i++);
                    if (childExp != null) {
                        val dx = ((x2+1)%2)*2-1;
                        val dy = ((y2+1)%2)*2-1;
                        val dz = ((z2+1)%2)*2-1;
                        this.multipoleExp.translateAndAddMultipole(scratch, scratch_array,
                          Vector3d(dx*halfSideLength, dy*halfSideLength, dz*halfSideLength),
                          myComplexK(dx,dy,dz), childExp, myWignerA((dx+1)/2, (dy+1)/2, (dz+1)/2));
                    }
                }
            }
        }

        sendMultipole(locallyEssentialTree, boxes, periodic);

        return this.multipoleExp;
    }

    protected def sendMultipole(locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)], periodic:Boolean) {
        // async send this box's multipole expansion to V-list
        if (vList != null) {
            val level = this.level;
            val x = this.x;
            val y = this.y;
            val z = this.z;
            val multipoleExp = this.multipoleExp;
            val thisLevelBoxes = boxes(level);
            val vListPlaces = new HashSet[Int]();
            for ([p] in vList) {
                vListPlaces.add(thisLevelBoxes.dist(vList(p)).id);
            }
            for(placeId in vListPlaces) {
                at(Place(placeId)) async {
                    val multipoleCopies = locallyEssentialTree().multipoleCopies(level);
                    if (periodic) {
                        val copyRegion = multipoleCopies.region;
                        val levelDim = Math.pow2(level) as Int;
                        for (xTrans in -1..1) {
                            for (yTrans in -1..1) {
                                for (zTrans in -1..1) {
                                    if (copyRegion.contains(x+xTrans*levelDim, y+yTrans*levelDim, z+zTrans*levelDim)) {
                                        multipoleCopies(x+xTrans*levelDim, y+yTrans*levelDim, z+zTrans*levelDim) = multipoleExp;
                                    }
                                }
                            }
                        }
                    } else {
                        multipoleCopies(x,y,z) = multipoleExp;
                    }
                }
            }
        }
    }

    protected def downward(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)], numLevels:Int, periodic:Boolean):Double {
        constructLocalExpansion(size, fmmOperators, parentLocalExpansion, locallyEssentialTree);

        val childLevel = level+1;
        val childLevelBoxes = boxes(childLevel);

        val parentExp = localExp;

        val farField = finish (SumReducer()) {
            for (i in 0..7) {
                val dx = i / 4;
                val dy = i % 4 / 2;
                val dz = i % 2;
                val childX = 2*this.x + dx;
                val childY = 2*this.y + dy;
                val childZ = 2*this.z + dz;

                at(childLevelBoxes.dist(childX,childY,childZ)) async {
                    val childBox = boxes(childLevel)(childX,childY,childZ);
                    if (childBox != null) {
                        offer childBox.downward(size, parentExp, fmmOperators, locallyEssentialTree, boxes, numLevels, periodic);
                    }
                }
            }
        };
        return farField;
    }

    protected def constructLocalExpansion(size:Double, fmmOperators:PlaceLocalHandle[FmmOperators], parentLocalExpansion:LocalExpansion, locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree]) {
        val myOperators = fmmOperators();
        val sideLength = size / Math.pow2(level);
        val myComplexK = myOperators.complexK;
        val myWignerB = myOperators.wignerB;
        val myWignerC = myOperators.wignerC;
        val thisLevelMultipoleCopies = locallyEssentialTree().multipoleCopies(level);

        // transform and add multipole expansions from same level
        localExp.terms.clear();
        val numTerms = localExp.p;
        val scratch = new MultipoleExpansion(numTerms);    
        val scratch_array = new Array[Complex](numTerms+1);
        val vList = getVList();
        for ([p] in vList) {
            
            val boxIndex2 = vList(p);
            // TODO - should be able to detect Point rank and inline
            val x2 = boxIndex2(0);
            val y2 = boxIndex2(1);
            val z2 = boxIndex2(2);
            val box2MultipoleExp = thisLevelMultipoleCopies(x2, y2, z2);
           
            if (box2MultipoleExp != null) {
               // TODO should be able to inline Point
                val dx2 = x2-this.x;
                val dy2 = y2-this.y;
                val dz2 = z2-this.z;
                localExp.transformAndAddToLocal(scratch, scratch_array,
			        Vector3d(dx2*sideLength, dy2*sideLength, dz2*sideLength), 
					myComplexK(dx2,dy2,dz2), box2MultipoleExp, myWignerB(dx2,dy2,dz2) );
            }
        }

        if (parentLocalExpansion != null) {
            // translate and add parent local expansion
            val dx = 2*(this.x%2)-1;
            val dy = 2*(this.y%2)-1;
            val dz = 2*(this.z%2)-1;

            localExp.translateAndAddLocal(scratch, scratch_array,
                Vector3d(dx*0.5*sideLength, dy*0.5*sideLength, dz*0.5*sideLength),
                myComplexK(dx,dy,dz), parentLocalExpansion, myWignerC((dx+1)/2, (dy+1)/2, (dz+1)/2));
        }

    }

    public def getVList() = this.vList;

    public def setVList(vList : Rail[Point(3)]) {
        this.vList = vList;
    }

    /**
     * Creates the V-list for this box.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent.
     */
    public def createVList(ws : Int) {
        val levelDim = Math.pow2(this.level);
        val xOffset = this.x%2 == 1 ? -1 : 0;
        val yOffset = this.y%2 == 1 ? -1 : 0;
        val zOffset = this.z%2 == 1 ? -1 : 0;
        val vList = new ArrayList[Point(3)]();
        for (x in Math.max(0,this.x-2*ws+xOffset)..Math.min(levelDim-1,this.x+2*ws+1+xOffset)) {
            for (y in Math.max(0,this.y-2*ws+yOffset)..Math.min(levelDim-1,this.y+2*ws+1+yOffset)) {
                for (z in Math.max(0,this.z-2*ws+zOffset)..Math.min(levelDim-1,this.z+2*ws+1+zOffset)) {
                    if (wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        this.vList = vList.toArray();
    }

    /**
     * Creates the V-list for this box for use with
     * the periodic FMM.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent.
     */
    public def createVListPeriodic(ws : Int) {
        val xOffset = this.x%2 == 1 ? -1 : 0;
        val yOffset = this.y%2 == 1 ? -1 : 0;
        val zOffset = this.z%2 == 1 ? -1 : 0;
        val vList = new ArrayList[Point(3)]();
        for (x in (this.x-2*ws+xOffset)..(this.x+2*ws+1+xOffset)) {
            for (y in (this.y-2*ws+yOffset)..(this.y+2*ws+1+yOffset)) {
                for (z in (this.z-2*ws+zOffset)..(this.z+2*ws+1+zOffset)) {
                    if (wellSeparated(ws, x, y, z)) {
                        vList.add(Point.make(x,y,z));
                    }
                }
            }
        }
        this.vList = vList.toArray();
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    public def toString(): String {
        return "FmmBox level " + level + " (" + x + "," + y + "," + z + ")";
    }
}

