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
import x10.util.HashSet;
import x10.util.Pair;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents a shared octant in the 3D division of space for the
 * Fast Multipole Method.
 * @author milthorpe
 */
public class SharedOctant extends Octant implements Comparable[SharedOctant] {
    public val children:Rail[Octant] = new Array[Octant](8);

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(id:OctantId, numTerms:Int, localData:FmmLocalData, dMax:UByte) {
        super(id, numTerms);
        val levelDim = (Math.pow2(dMax) / Math.pow2(id.level));
        var i:Int = 0;
        for (x2 in (2*id.x)..(2*id.x+1)) {
            for (y2 in (2*id.y)..(2*id.y+1)) {
                for (z2 in (2*id.z)..(2*id.z+1)) {
                    val childOctantId = OctantId(x2 as UByte, y2 as UByte, z2 as UByte, id.level+1UY);
                    val placeId = localData.getPlaceId(childOctantId.getAnchor(dMax));
                    if (placeId != here.id && placeId >= 0 && placeId < Place.MAX_PLACES) {
                        // the child octant is not held at this place
                        //Console.OUT.println("at " + here + " octant " + id + " creating ghost for octant " + childOctantId + " held at " + placeId);
                        children(i) = new GhostOctant(childOctantId, placeId);
                    }
                    i++;
                 }
            }
        }
    }

    public def compareTo(b:SharedOctant):Int = id.compareTo(b.id);

    /** 
     * For each shared octant, combines multipole expansions for <= 8 child
     * octants into a single multipole expansion for the parent octant.
     */
    protected def upward(localData:PlaceLocalHandle[FmmLocalData], size:Double, dMax:UByte, periodic:Boolean):Pair[Int,MultipoleExpansion] {
        //Console.OUT.println("at " + here + " SharedOctant.upward for " + id + " children.size = " + children.size);
        numAtoms = 0; // reset
        this.parentAdded = false;

        val childExpansions = new Array[Pair[Int,MultipoleExpansion]](8);
        finish {
            for (i in children) {
                val child = children(i);
                if (child != null) async {
                    childExpansions(i) = child.upward(localData, size, dMax, periodic);
                }
            }
            multipoleExp.terms.clear();
            localExp.terms.clear();
        }

        val fmmOperators = localData().fmmOperators;
        val myComplexK = fmmOperators.complexK;
        val myWignerA = fmmOperators.wignerA;
        val halfSideLength = size / Math.pow2(id.level+2);
        val numTerms = multipoleExp.p;
        val scratch = new MultipoleExpansion(numTerms);    
        val scratch_array = new Array[Complex](numTerms+1);
        var i:Int=0;
        var nonNullChildren:Boolean = false;
        for (x2 in (2*id.x)..(2*id.x+1)) {
            for (y2 in (2*id.y)..(2*id.y+1)) {
                for (z2 in (2*id.z)..(2*id.z+1)) {
                    val childExp = childExpansions(i++);
                    numAtoms = numAtoms + childExp.first;
                    if (childExp.second != null) {
                        nonNullChildren = true;
                        val dx = ((x2+1)%2)*2-1;
                        val dy = ((y2+1)%2)*2-1;
                        val dz = ((z2+1)%2)*2-1;
                        this.multipoleExp.translateAndAddMultipole(scratch, scratch_array,
                          Vector3d(dx*halfSideLength, dy*halfSideLength, dz*halfSideLength),
                          myComplexK(dx,dy,dz), childExp.second, myWignerA((dx+1)/2, (dy+1)/2, (dz+1)/2));
                    }
                }
            }
        }

        atomic this.multipoleReady = true;

        if (nonNullChildren) {
            sendMultipole(localData, dMax, periodic);
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
        return Pair[Int,MultipoleExpansion](numAtoms, this.multipoleExp);
    }

    protected def downward(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, dMax:UByte, periodic:Boolean):Double {
        //Console.OUT.println("at " + here + " SharedOctant.downward for " + id + " numAtoms = " + numAtoms);

        this.multipoleReady = false; // reset

        if (numAtoms > 0) {
            constructLocalExpansion(localData, size, parentLocalExpansion);

            val parentExp = localExp;

            val farField = finish (SumReducer()) {
                for (i in 0..(children.size-1)) {
                    val childOctant = children(i);
                    if (childOctant != null) async {
                        val dx = i / 4;
                        val dy = i % 4 / 2;
                        val dz = i % 2;
                        val childX = 2*id.x + dx;
                        val childY = 2*id.y + dy;
                        val childZ = 2*id.z + dz;
                        offer childOctant.downward(localData, size, parentExp, dMax, periodic);
                    }
                }
            };

            return farField;

        } else {
            return 0.0;
        }
    }

    public def addToCombinedVSet(combinedVSet:HashSet[OctantId], ws:Int) {
        super.addToCombinedVSet(combinedVSet, ws);
        for (i in 0..(children.size-1)) {
            val childOctant = children(i);
            if (children(i) != null) {
                childOctant.addToCombinedVSet(combinedVSet, ws);
            }
        }
    }

    public def getDescendant(octantId:OctantId):Octant {
        //Console.OUT.println("at " + here + " SharedOctant.getDescendant(" + octantId + ") on " + this.id);
        if (octantId == this.id) return this;
        for (i in 0..(children.size-1)) {
            val childOctant = children(i);
            if (childOctant != null) {
                val desc = childOctant.getDescendant(octantId);
                if (desc != null) {
                    return desc;
                }
            }
        }
        return null;
    }

    public def toString(): String {
        return "SharedOctant " + id;
    }
}

