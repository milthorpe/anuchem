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
import x10.util.ArrayUtils;
import x10.util.HashSet;
import x10.util.Pair;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents an octant in the 3D division of space
 * for the fast multipole method.
 * @author milthorpe
 */
public abstract class Octant implements Comparable[Octant] {
    public id:OctantId;

    public var parent:Octant;

    /** The number of atoms in all boxes below this box. */
    public var numAtoms:Int;

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private var vList:Rail[OctantId];

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp:MultipoleExpansion;

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp:LocalExpansion;

    /** 
     * Flag set to true when octant multipole expansion is consistent i.e.
     * completed in upward pass, and set to false again during downward pass.
     */
    public var multipoleReady:Boolean = false;

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(id:OctantId, numTerms:Int) {
        this.id = id;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
    }

    public def compareTo(b:Octant):Int = id.compareTo(b.id);

    public def getCentre(size:Double):Point3d {
        dim:Int = Math.pow2(id.level);
        sideLength:Double = size / dim;
        offset:Double = 0.5 * size;
        return Point3d( (id.x + 0.5) * sideLength - offset,
                        (id.y + 0.5) * sideLength - offset,
                        (id.z + 0.5) * sideLength - offset);
    }

    abstract protected def downward(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, dMax:UByte):Double;

    /** 
     * Generates and combines multipole expansions for all descendants into an
     * expansion for this box. Note: non-blocking - the top-level call to this 
     * method must be enclosed in a finish statement.
     */
    abstract protected def upward(localData:PlaceLocalHandle[FmmLocalData], size:Double, dMax:UByte):Pair[Int,MultipoleExpansion];

    protected def constructLocalExpansion(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion) {
        val local = localData();
        val sideLength = size / Math.pow2(id.level);
        val myComplexK = local.fmmOperators.complexK;
        val myWignerB = local.fmmOperators.wignerB;
        val myWignerC = local.fmmOperators.wignerC;
        val locallyEssentialTree = local.locallyEssentialTree;

        // transform and add multipole expansions from same level
        val numTerms = localExp.p;
        val scratch = new MultipoleExpansion(numTerms);    
        val scratch_array = new Array[Complex](numTerms+1);
        val vList = getVList();
        for ([p] in vList) {
            val octantIndex2 = vList(p);
            val box2MultipoleExp = locallyEssentialTree.getMultipoleForOctant(octantIndex2);
           
            if (box2MultipoleExp != null) {
                //Console.OUT.println("at " + here + " adding multipole for " + octantIndex2 + " to " + id);
                val dx2 = (octantIndex2.x as Int)-id.x;
                val dy2 = (octantIndex2.y as Int)-id.y;
                val dz2 = (octantIndex2.z as Int)-id.z;
                localExp.transformAndAddToLocal(scratch, scratch_array,
			        Vector3d(dx2*sideLength, dy2*sideLength, dz2*sideLength), 
					myComplexK(dx2,dy2,dz2), box2MultipoleExp, myWignerB(dx2,dy2,dz2) );
                //Console.OUT.println("at " + here + " added multipole for " + octantIndex2 + " to " + id);
            }
        }

        if (parentLocalExpansion != null) {
            // translate and add parent local expansion
            val dx = 2*(id.x%2)-1;
            val dy = 2*(id.y%2)-1;
            val dz = 2*(id.z%2)-1;

            localExp.translateAndAddLocal(scratch, scratch_array,
                Vector3d(dx*0.5*sideLength, dy*0.5*sideLength, dz*0.5*sideLength),
                myComplexK(dx,dy,dz), parentLocalExpansion, myWignerC((dx+1)/2, (dy+1)/2, (dz+1)/2));
        }
    }

    protected def sendMultipole(localData:PlaceLocalHandle[FmmLocalData], dMax:UByte) {
        // async send this box's multipole expansion to V-list
        //Console.OUT.println("at " + here + " sending multipole for " + id);
        if (vList != null) {
            val local = localData();
            val id = this.id;
            val multipoleExp = this.multipoleExp;
            val vListPlaces = new HashSet[Int]();
            for ([p] in vList) {
                val placeId = local.getPlaceId(vList(p).getAnchor(dMax));
                if (placeId >= 0 && placeId < Place.MAX_PLACES) {
                    //Console.OUT.println("at " + here + " sending multipole for " + id + " " + vList(p) + " held at " + placeId);
                    vListPlaces.add(placeId);
                }
            }
            for(placeId in vListPlaces) {
                at(Place(placeId)) async {
                    //Console.OUT.println("at " + here + " sending multipole for " + id + " to place " + placeId);
                    localData().locallyEssentialTree.setMultipoleForOctant(id, multipoleExp);
                }
            }
        }
    }

    /**
     * Creates the V-list for this box.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent.
     */
    public def createVList(ws:Int) {
        val levelDim = Math.pow2(id.level);
        val xOffset = id.x%2 == 1UY ? -1 : 0;
        val yOffset = id.y%2 == 1UY ? -1 : 0;
        val zOffset = id.z%2 == 1UY ? -1 : 0;
        val vList = new ArrayList[OctantId]();
        for (x in Math.max(0,id.x-2*ws+xOffset)..Math.min(levelDim-1,id.x+2*ws+1+xOffset)) {
            val x2 = x as UByte;
            for (y in Math.max(0,id.y-2*ws+yOffset)..Math.min(levelDim-1,id.y+2*ws+1+yOffset)) {
                val y2 = y as UByte;
                for (z in Math.max(0,id.z-2*ws+zOffset)..Math.min(levelDim-1,id.z+2*ws+1+zOffset)) {
                    val z2 = z as UByte;
                    if (wellSeparated(ws, x, y, z)) {
                        vList.add(OctantId(x2,y2,z2,id.level));
                    }
                }
            }
        }
        this.vList = vList.toArray();
    }

    public def getVList() = this.vList;

    public def addToCombinedVSet(combinedVSet:HashSet[OctantId], ws:Int) {
        createVList(ws);
        for ([p] in vList) {
            combinedVSet.add(vList(p));
        }
    }

    /**
     * Returns true if this box is well-separated from <code>x,y,z</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    public def wellSeparated(ws:Int, x2:Int, y2:Int, z2:Int) : Boolean {
        return Math.abs(id.x - x2) > ws 
            || Math.abs(id.y - y2) > ws 
            || Math.abs(id.z - z2) > ws;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }

    public def toString(): String {
        return "Octant " + id;
    }
}

