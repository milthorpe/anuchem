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
import x10.util.HashSet;
import x10.util.Pair;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents a parent octant in the 3D division of space for the
 * Fast Multipole Method.
 * @author milthorpe
 */
public class ParentOctant extends Octant implements Comparable[ParentOctant] {

    public val children:Rail[Octant] = new Array[Octant](8);

    /** The number of atoms in all boxes below this box. */
    private var numAtoms:Int;

    /**
     * Creates a new FmmBox with multipole and local expansions
     * of the given number of terms.
     */
    public def this(id:OctantId, numTerms:Int, ws:Int, dMax:UByte) {
        super(id, numTerms, ws, dMax);
    }

    public def countOctants():Int {
        var octants:Int = 0;
        for ([i] in children) {
            val child = children(i);
            if (child != null) {
                octants += child.countOctants();
            }
        }
        return octants;
    }

    public def ghostOctants():Int {
        var ghostOctants:Int = 0;
        for ([i] in children) {
            val child = children(i);
            if (child != null) {
                ghostOctants += child.ghostOctants();
            }
        }
        return ghostOctants;
    }

    public def compareTo(b:ParentOctant):Int = id.compareTo(b.id);

    public def numAtoms() = numAtoms;

    /** 
     * For each shared octant, combines multipole expansions for <= 8 child
     * octants into a single multipole expansion for the parent octant.
     */
    protected def upward():Pair[Int,MultipoleExpansion] {
        //Console.OUT.println("at " + here + " ParentOctant.upward for " + id + " children.size = " + children.size);
        numAtoms = 0; // reset

        val childExpansions = new Array[Pair[Int,MultipoleExpansion]](8);
        finish {
            for ([i] in children) {
                val child = children(i);
                if (child != null) async {
                    childExpansions(i) = child.upward();
                }
            }
            multipoleExp.clear();
            localExp.clear();
        }

        val local = FastMultipoleMethod.localData;
        val fmmOperators = local.fmmOperators;
        val myComplexK = fmmOperators.complexK;
        val myWignerA = fmmOperators.wignerA;
        val halfSideLength = local.size / Math.pow2(id.level+2);
        val numTerms = multipoleExp.p;
        val scratch = FmmScratch.getWorkerLocal();
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
                        this.multipoleExp.translateAndAddMultipole(scratch.exp, scratch.array,
                          Vector3d(dx*halfSideLength, dy*halfSideLength, dz*halfSideLength),
                          myComplexK(dx,dy,dz), childExp.second, myWignerA((dx+1)/2, (dy+1)/2, (dz+1)/2));
                    }
                }
            }
        }

        atomic this.multipoleReady = true;

        if (nonNullChildren) {
            sendMultipole();
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
        return Pair[Int,MultipoleExpansion](numAtoms, this.multipoleExp);
    }

    public def multipolesToLocal() {
        for ([i] in children) {
            val child = children(i);
            if (child != null && !(child instanceof GhostOctant)) {
                async child.multipolesToLocal();
            }
        }
        super.multipolesToLocal();
    }

    protected def downward(parentLocalExpansion:LocalExpansion):Double {
        //Console.OUT.println("at " + here + " ParentOctant.downward for " + id + " numAtoms = " + numAtoms);

        this.multipoleReady = false; // reset

        if (numAtoms > 0) {
            addParentExpansion(parentLocalExpansion);

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
                        offer childOctant.downward(parentExp);
                    }
                }
            };
            return farField;
        } else {
            return 0.0;
        }
    }

    public def toString(): String {
        return "ParentOctant " + id;
    }
}

