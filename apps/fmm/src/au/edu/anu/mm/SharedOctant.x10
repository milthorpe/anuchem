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
    public def this(id:OctantId, numTerms:Int) {
        super(id, numTerms);
    }

    public def compareTo(b:SharedOctant):Int = id.compareTo(b.id);

    /** 
     * For each shared octant, combines multipole expansions for <= 8 child
     * octants into a single multipole expansion for the parent octant.
     */
    protected def upward(localData:PlaceLocalHandle[FmmLocalData], size:Double, periodic:Boolean):Pair[Int,MultipoleExpansion] {
        Console.OUT.println("SharedOctant.upward for " + id);
        numAtoms = 0; // reset

        val childExpansions = new Array[Pair[Int,MultipoleExpansion]](8);
        finish {
            for (i in 0..(children.size-1)) {
                val childOctant = children(i);
                if (childOctant != null) {
                    async childExpansions(i) = childOctant.upward(localData, size, periodic);
                }
                multipoleExp.terms.clear();
            }
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

        if (nonNullChildren) {
            sendMultipole(localData, periodic);
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
        return Pair[Int,MultipoleExpansion](numAtoms, this.multipoleExp);
    }

    protected def downward(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, numLevels:Int, periodic:Boolean):Double {
        Console.OUT.println("SharedOctant.downward for " + id + " numAtoms = " + numAtoms);
        if (numAtoms > 0) {
            constructLocalExpansion(localData, size, parentLocalExpansion);

            val parentExp = localExp;

            val farField = finish (SumReducer()) {
                for (i in 0..(children.size-1)) {
                    val dx = i / 4;
                    val dy = i % 4 / 2;
                    val dz = i % 2;
                    val childX = 2*id.x + dx;
                    val childY = 2*id.y + dy;
                    val childZ = 2*id.z + dz;

                    val childOctant = children(i);
                    if (childOctant != null) {
                        offer childOctant.downward(localData, size, parentExp, numLevels, periodic);
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

    public def getDescendant(id:OctantId):Octant {
        for (i in 0..(children.size-1)) {
            val childOctant = children(i);
            if (childOctant != null) {
                val desc = childOctant.getDescendant(id);
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

