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

import x10.util.ArrayList;
import x10.util.RailUtils;
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

    /** The multipole expansion of the charges within this box. */
    public val multipoleExp:MultipoleExpansion;

    /** The Taylor expansion of the potential within this box due to particles in well separated boxes. */
    public val localExp:LocalExpansion;

    /** The V-list for this octant, or null if this is a ghost octant. */
    public var vList:VList;

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
    }

    /**
     * Creates a new Octant with no expansions (used by GhostOctant).
     */
    protected def this(id:OctantId) {
        this.id = id;
        this.multipoleExp = null;
        this.localExp = null;
    }

    public def compareTo(b:Octant):Int = id.compareTo(b.id);

    public abstract def numAtoms():Long;
    
    public def getCentre(size:Double):Point3d {
        val dim = 1UY << id.level;
        val sideLength = size / dim;
        val offset = 0.5 * size;
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
    abstract protected def upward():Pair[Long,MultipoleExpansion];

    /*
     * Transform and add multipole expansions from same level
     */
    public def multipolesToLocal() {
        val local = FastMultipoleMethod.localData;
        val sideLength = local.size / (1UY << id.level);
        val myComplexK = local.fmmOperators.complexK;
        val myWignerB = local.fmmOperators.wignerB;
        val locallyEssentialTree = local.locallyEssentialTree;


        val numTerms = localExp.p;
        val scratch = FmmScratch.getWorkerLocal();
        for (octantIndex2 in vList) {
            val box2MultipoleExp = locallyEssentialTree.getMultipoleForOctant(octantIndex2.getMortonId());
           
            if (box2MultipoleExp != null) {
                //Console.OUT.println("at " + here + " adding multipole for " + octantIndex2 + " to " + id);
                val dx2 = (octantIndex2.x as Long)-id.x;
                val dy2 = (octantIndex2.y as Long)-id.y;
                val dz2 = (octantIndex2.z as Long)-id.z;
                localExp.transformAndAddToLocal(scratch.exp, scratch.array,
			        Vector3d(dx2*sideLength, dy2*sideLength, dz2*sideLength), 
					myComplexK(dx2,dy2,dz2), box2MultipoleExp, myWignerB(dx2,dy2,dz2) );
            }
        }
    }

    protected def addParentExpansion(parentLocalExpansion:LocalExpansion) {
        if (parentLocalExpansion != null) {
            val local = FastMultipoleMethod.localData;
            val sideLength = local.size / (1UY << id.level);
            val myComplexK = local.fmmOperators.complexK;
            val myWignerC = local.fmmOperators.wignerC;

            val numTerms = localExp.p;
            val scratch = FmmScratch.getWorkerLocal();
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
        if (vList != null) {
        // async send this box's multipole expansion to V-list
        //Console.OUT.println("at " + here + " sending multipole for " + id);
            val local = FastMultipoleMethod.localData;
            val mortonId = id.getMortonId();
            val multipoleExp = this.multipoleExp;
            val vListPlaces = new HashSet[Long]();
            for (octantId in vList) {
                val placeId = local.getPlaceId(octantId.getAnchor(local.dMax));
                if (placeId >= 0 && placeId < Place.numPlaces()) {
                    //Console.OUT.println("at " + here + " sending multipole for " + id + " to " + octantId + " held at " + placeId);
                    vListPlaces.add(placeId);
                }
            }
            for(placeId in vListPlaces) {
                if (placeId == here.id) {
                    local.locallyEssentialTree.setMultipoleForOctant(mortonId, multipoleExp);
                } else {
                    at(Place(placeId)) async {
                        FastMultipoleMethod.localData.locallyEssentialTree.setMultipoleForOctant(mortonId, multipoleExp);
                    }
                }
            }
        }
    }

    /** Estimates the number of octants in an octant's V-list, given the octant id. */
    public static def estimateVListSize(mortonId:UInt):Int {
        val id = OctantId.getFromMortonId(mortonId);
        return estimateVListSize(id, 1n);
    }

    /** Estimates the number of octants in an octant's V-list, given the octant id. */
    public static def estimateVListSize(id:OctantId, ws:Int):Int {
        val near = 2n*ws+1n;
        val farHalo = 4n*ws+2n;
        val maxExtent = (1n << id.level) - 1n;
        val x = id.x as Int;
        val y = id.y as Int;
        val z = id.z as Int;
        val xOffset = x%2n == 1n ? -1n : 0n;
        val yOffset = y%2n == 1n ? -1n : 0n;
        val zOffset = z%2n == 1n ? -1n : 0n;
        val nearMinX = Math.max(0n, x-ws);
        val nearMaxX = Math.min(maxExtent, x+ws);
        val nearMinY = Math.max(0n, y-ws);
        val nearMaxY = Math.min(maxExtent, y+ws);
        val nearMinZ = Math.max(0n, z-ws);
        val nearMaxZ = Math.min(maxExtent, z+ws);
        
        val neighbours = (nearMaxX-nearMinX+1n) * (nearMaxY-nearMinY+1n) * (nearMaxZ-nearMinZ+1n);

        val farExtent = 2n*ws;

        val farMinX = Math.max(0n, x+xOffset-farExtent);
        val farMaxX = Math.min(maxExtent, x+xOffset+farExtent+1n);
        val farMinY = Math.max(0n, y+yOffset-farExtent);
        val farMaxY = Math.min(maxExtent, y+yOffset+farExtent+1n);
        val farMinZ = Math.max(0n, z+zOffset-farExtent);
        val farMaxZ = Math.min(maxExtent, z+zOffset+farExtent+1n);
        
        val colleagues = (farMaxX-farMinX+1n) * (farMaxY-farMinY+1n) * (farMaxZ-farMinZ+1n);

        return colleagues - neighbours;
    }

    public def createVList(ws:Int, dMax:UByte) {
        vList = new VList(id, ws, dMax);
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

    public def toString(): String {
        return "Octant " + id;
    }

    /** 
     * The V-list consists of the children of those boxes 
     * not well-separated from this box's parent.
     */
    private class VList implements Iterable[OctantId] {
        val ws:Int;
        val level:UByte;
        val minX:UByte;
        val maxX:UByte;
        val minY:UByte;
        val maxY:UByte;
        val minZ:UByte;
        val maxZ:UByte;

        public def this(id:OctantId, ws:Int, dMax:UByte) {
            val levelDim = 1UY << id.level;
            val xOffset = id.x%2UY == 1UY ? -1 : 0;
            val yOffset = id.y%2UY == 1UY ? -1 : 0;
            val zOffset = id.z%2UY == 1UY ? -1 : 0;
            val extent = 2*ws;
            minX = Math.max(0,id.x+xOffset-extent) as UByte;
            maxX = Math.min(levelDim-1,id.x+xOffset+extent+1) as UByte;
            minY = Math.max(0,id.y+yOffset-extent) as UByte;
            maxY = Math.min(levelDim-1,id.y+yOffset+extent+1) as UByte;
            minZ = Math.max(0,id.z+zOffset-extent) as UByte;
            maxZ = Math.min(levelDim-1,id.z+zOffset+extent+1) as UByte;
            this.ws = ws;
            this.level = id.level;
        }

        /**
         * Returns true if this box is well-separated from <code>x,y,z</code>
         * on the same level, i.e. if there are at least <code>ws</code>
         * boxes separating them.
         */
        private def wellSeparated(ws:Int, x2:Int, y2:Int, z2:Int) : Boolean {
            return Math.abs(Octant.this.id.x - x2) > ws 
                || Math.abs(Octant.this.id.y - y2) > ws 
                || Math.abs(Octant.this.id.z - z2) > ws;
        }

        public def iterator() = new VListIterator();

        public class VListIterator implements Iterator[OctantId] {
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
                    if (wellSeparated(ws, x, y, z)) {
                        return true;
                    } else {
                        moveToNext();
                        if (x <= maxX && wellSeparated(ws, x, y, z)) return true;
                    }
                }
                return false;
            }
            
            public def next():OctantId {
                if (x <= maxX && wellSeparated(ws, x, y, z)) {
                    val res = new OctantId(x, y, z, level);
                    moveToNext();
                    return res;
                } else {
                    throw new UnsupportedOperationException("reached end of vList for " + Octant.this.id);
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
                } while(x <= maxX && !wellSeparated(ws, x, y, z));
            }
        }
    }
}

