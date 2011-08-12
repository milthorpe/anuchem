/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.util.ArrayList;

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

    protected def downward(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:PlaceLocalHandle[FmmOperators], locallyEssentialTree:PlaceLocalHandle[LocallyEssentialTree], boxes:Rail[DistArray[FmmBox](3)]):Double {
        val myOperators = fmmOperators();
        addMultipolesAtSameLevel(size, myOperators, locallyEssentialTree());
        addParentExpansion(size, parentLocalExpansion, myOperators);

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

                async at(childLevelBoxes.dist(childX,childY,childZ)) {
                    val childBox = boxes(childLevel)(childX,childY,childZ);
                    if (childBox != null) {
                        offer childBox.downward(size, parentExp, fmmOperators, locallyEssentialTree, boxes);
                    }
                }
            }
        };
        return farField;
    }

    protected def addMultipolesAtSameLevel(size:Double, fmmOperators:FmmOperators, locallyEssentialTree:LocallyEssentialTree) {
        val sideLength = size / Math.pow2(level);
        val myComplexK = fmmOperators.complexK;
        val myWignerB = fmmOperators.wignerB;
        val thisLevelMultipoleCopies = locallyEssentialTree.multipoleCopies(level);

        // transform and add multipole expansions from same level
        val numTerms = localExp.p;
        val scratch = new MultipoleExpansion(numTerms);    
        val scratch_array = new Array[Complex](-numTerms..numTerms) as Array[Complex](1){rect,rail==false};
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
    }

    protected def addParentExpansion(size:Double, parentLocalExpansion:LocalExpansion, fmmOperators:FmmOperators) {
        val sideLength = size / Math.pow2(level);
        val myComplexK = fmmOperators.complexK;
        val myWignerC = fmmOperators.wignerC;

        if (parentLocalExpansion != null) {
            // translate and add parent local expansion
            val dx = 2*(this.x%2)-1;
            val dy = 2*(this.y%2)-1;
            val dz = 2*(this.z%2)-1;

            localExp.translateAndAddLocal(
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

