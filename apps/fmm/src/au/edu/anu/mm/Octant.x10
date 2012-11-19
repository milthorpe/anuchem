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

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private val vList:Rail[OctantId];

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
     * Creates a new Octant with multipole and local expansions
     * of the given number of terms.
     */
    public def this(id:OctantId, numTerms:Int, ws:Int, dMax:UByte) {
        this.id = id;
        this.multipoleExp = new MultipoleExpansion(numTerms);
        this.localExp = new LocalExpansion(numTerms);
        this.vList = createVList(ws, dMax);
    }

    /**
     * Creates a new Octant with no expansions (used by GhostOctant).
     */
    protected def this(id:OctantId) {
        this.id = id;
        this.multipoleExp = null;
        this.localExp = null;
        this.vList = null;
    }

    public def compareTo(b:Octant):Int = id.compareTo(b.id);

    public abstract def numAtoms():Int;
    
    public def getCentre(size:Double):Point3d {
        dim:Int = Math.pow2(id.level);
        sideLength:Double = size / dim;
        offset:Double = 0.5 * size;
        return Point3d( (id.x + 0.5) * sideLength - offset,
                        (id.y + 0.5) * sideLength - offset,
                        (id.z + 0.5) * sideLength - offset);
    }

    abstract public def countOctants():Int;
    abstract public def ghostOctants():Int;

    abstract protected def downward(parentLocalExpansion:LocalExpansion):Double;

    /** 
     * Generates and combines multipole expansions for all descendants into an
     * expansion for this box. Note: non-blocking - the top-level call to this 
     * method must be enclosed in a finish statement.
     */
    abstract protected def upward():Pair[Int,MultipoleExpansion];

    protected def constructLocalExpansion(parentLocalExpansion:LocalExpansion) {
        val local = FastMultipoleMethod.localData;
        val sideLength = local.size / Math.pow2(id.level);
        val myComplexK = local.fmmOperators.complexK;
        val myWignerB = local.fmmOperators.wignerB;
        val myWignerC = local.fmmOperators.wignerC;
        val locallyEssentialTree = local.locallyEssentialTree;

        // transform and add multipole expansions from same level
        val numTerms = localExp.p;
        val scratch = FmmScratch.getWorkerLocal();
        val vList = getVList();
        //Console.OUT.println(id + " vList " + vList.size);
        for ([p] in vList) {
            val octantIndex2 = vList(p);
            val box2MultipoleExp = locallyEssentialTree.getMultipoleForOctant(octantIndex2.getMortonId());
           
            if (box2MultipoleExp != null) {
                //Console.OUT.println("at " + here + " adding multipole for " + octantIndex2 + " to " + id);
                val dx2 = (octantIndex2.x as Int)-id.x;
                val dy2 = (octantIndex2.y as Int)-id.y;
                val dz2 = (octantIndex2.z as Int)-id.z;
                localExp.transformAndAddToLocal(scratch.exp, scratch.array,
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

            localExp.translateAndAddLocal(scratch.exp, scratch.array,
                Vector3d(dx*0.5*sideLength, dy*0.5*sideLength, dz*0.5*sideLength),
                myComplexK(dx,dy,dz), parentLocalExpansion, myWignerC((dx+1)/2, (dy+1)/2, (dz+1)/2));
        }
    }

    protected def sendMultipole() {
        // async send this box's multipole expansion to V-list
        //Console.OUT.println("at " + here + " sending multipole for " + id);
        if (vList != null) {
            val local = FastMultipoleMethod.localData;
            val mortonId = id.getMortonId();
            val multipoleExp = this.multipoleExp;
            val vListPlaces = new HashSet[Int]();
            for ([p] in vList) {
                val placeId = local.getPlaceId(vList(p).getAnchor(local.dMax));
                if (placeId >= 0 && placeId < Place.MAX_PLACES) {
                    //Console.OUT.println("at " + here + " sending multipole for " + id + " " + vList(p) + " held at " + placeId);
                    vListPlaces.add(placeId);
                }
            }
            for(placeId in vListPlaces) {
                if (placeId == here.id) {
                    local.locallyEssentialTree.setMultipoleForOctant(mortonId, multipoleExp);
                } else {
                    at(Place(placeId)) async {
                        //Console.OUT.println("at " + here + " sending multipole for " + id + " to place " + placeId);
                        FastMultipoleMethod.localData.locallyEssentialTree.setMultipoleForOctant(mortonId, multipoleExp);
                    }
                }
            }
        }
    }

    /** Estimates the number of octants in an octant's V-list, given the octant id. */
    public static def estimateVListSize(mortonId:UInt):Int {
        val id = OctantId.getFromMortonId(mortonId);
        return estimateVListSize(id, 1);
    }

    /** Estimates the number of octants in an octant's V-list, given the octant id. */
    public static def estimateVListSize(id:OctantId, ws:Int):Int {
        val near = 2*ws+1;
        val farHalo = 4*ws+2;
        val maxExtent = (1 << id.level) - 1;
        val x = id.x as Int;
        val y = id.y as Int;
        val z = id.z as Int;
        val xOffset = x%2 == 1 ? -1 : 0;
        val yOffset = y%2 == 1 ? -1 : 0;
        val zOffset = z%2 == 1 ? -1 : 0;
        val nearMinX = Math.max(0, x-ws);
        val nearMaxX = Math.min(maxExtent, x+ws);
        val nearMinY = Math.max(0, y-ws);
        val nearMaxY = Math.min(maxExtent, y+ws);
        val nearMinZ = Math.max(0, z-ws);
        val nearMaxZ = Math.min(maxExtent, z+ws);
        
        val neighbours = (nearMaxX-nearMinX+1) * (nearMaxY-nearMinY+1) * (nearMaxZ-nearMinZ+1);

        val farExtent = 2*ws;

        val farMinX = Math.max(0, x+xOffset-farExtent);
        val farMaxX = Math.min(maxExtent, x+xOffset+farExtent+1);
        val farMinY = Math.max(0, y+yOffset-farExtent);
        val farMaxY = Math.min(maxExtent, y+yOffset+farExtent+1);
        val farMinZ = Math.max(0, z+zOffset-farExtent);
        val farMaxZ = Math.min(maxExtent, z+zOffset+farExtent+1);
        
        val colleagues = (farMaxX-farMinX+1) * (farMaxY-farMinY+1) * (farMaxZ-farMinZ+1);

        return colleagues - neighbours;
    }

    /**
     * Returns a cost estimate per interaction (in ns) of V-list calculation.
     */
    public static def estimateVListCost(numTerms:Int):Long {
        val local = FastMultipoleMethod.localData;
        val myComplexK = local.fmmOperators.complexK;
        val myWignerB = local.fmmOperators.wignerB;

        val scratch = FmmScratch.getWorkerLocal();
        val localExp = new LocalExpansion(numTerms);
        val randomExp = new MultipoleExpansion(numTerms);
        val start = System.nanoTime();
        for (i in 1..10) {
            localExp.transformAndAddToLocal(scratch.exp, scratch.array,
		        Vector3d(1.0, 1.0, 1.0), 
				myComplexK(1,1,1), randomExp, myWignerB(1,1,1) );
        }
        val stop = System.nanoTime();
        localExp.clear();
        return (stop-start)/10L;
    }

    /**
     * Creates the V-list for this box.
     * The V-list consists of the children of those boxes not 
     * well-separated from the parent.
     */
    private def createVList(ws:Int, dMax:UByte) {
        val vList = new Array[OctantId](estimateVListSize(id, ws));
        val levelDim = Math.pow2(id.level);
        val xOffset = id.x%2 == 1UY ? -1 : 0;
        val yOffset = id.y%2 == 1UY ? -1 : 0;
        val zOffset = id.z%2 == 1UY ? -1 : 0;
        val extent = 2*ws;
        var i:Int=0;
        for (x in Math.max(0,id.x+xOffset-extent)..Math.min(levelDim-1,id.x+xOffset+extent+1)) {
            val x2 = x as UByte;
            for (y in Math.max(0,id.y+yOffset-extent)..Math.min(levelDim-1,id.y+yOffset+extent+1)) {
                val y2 = y as UByte;
                for (z in Math.max(0,id.z+zOffset-extent)..Math.min(levelDim-1,id.z+zOffset+extent+1)) {
                    val z2 = z as UByte;
                    if (wellSeparated(ws, x, y, z)) {
                        vList(i++) = OctantId(x2,y2,z2,id.level);
                    }
                }
            }
        }
        return vList;
    }

    public def getVList() = this.vList;

    /**
     * Returns true if this box is well-separated from <code>x,y,z</code>
     * on the same level, i.e. if there are at least <code>ws</code>
     * boxes separating them.
     */
    private def wellSeparated(ws:Int, x2:Int, y2:Int, z2:Int) : Boolean {
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

