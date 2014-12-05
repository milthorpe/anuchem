/*
 *  This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.chem.mm;

import x10.compiler.Inline;
import x10.util.ArrayList;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * Holds all particle data for a given place in a molecular mechanics
 * simulation.  Atom features (atomType, position, velocity, force) are held
 * in separate Rails, indexed by a 'local atom index' at this place.
 * The globalIndex for an atom is unique across all places.
 */
public class ParticleData {
    public var description:String;

    public var atomTypes:Rail[AtomType];

    /** Global index of each atom */
    public val globalIndex = new ArrayList[Long]();
    /** Atom type index for each atom */
    public val atomTypeIndex = new ArrayList[Int]();
    /** Atom centers in nm */
    public val x = new ArrayList[Point3d]();
    /** Atom velocities */
    public val dx = new ArrayList[Vector3d]();
    /** Instantaneous force on atoms */
    public val fx = new ArrayList[Vector3d]();

    public final def numAtoms() = globalIndex.size();

    public def addAtom(index:Long, atomType:Int, center:Point3d, velocity:Vector3d) {
        globalIndex.add(index);
        atomTypeIndex.add(atomType);
        x.add(center);
        dx.add(velocity);
        fx.add(Vector3d.NULL);
    }

    public def addAtom(index:Long, atomType:Int, center:Point3d) {
        addAtom(index, atomType, center, Vector3d.NULL);
    }

    public def removeAtom(index:Long) {
        globalIndex.removeAt(index);
        atomTypeIndex.removeAt(index);
        x.removeAt(index);
        dx.removeAt(index);
        fx.removeAt(index);
    }

    /**
     * Sorts a portion of the atoms into ascending order.
     * @param lo the index of the start of the range to be sorted
     * @param hi the index of the end of the range to be sorted
     * @param cmp the comparison function to use
     */
    public def sortAtoms(lo:Long, hi:Long, cmp:(Long,Long)=>Int) {
        if (hi <= lo) return;
        var l:Long = lo - 1;
        var h:Long = hi;
        while (true) {
            while (cmp(++l,hi)<0);
            while (cmp(hi, --h)<0 && h>lo);
            if (l >= h) break;
            swapAtoms(l, h);
        }
        swapAtoms(l, hi);
        sortAtoms(lo, l-1, cmp);
        sortAtoms(l+1, hi, cmp);
    }

    /**
     * Sorts a portion of the atoms into ascending order based on an external
     * feature array.
     * @param a the Rail containing the feature upon which to sort
     * @param lo the index of the start of the range to be sorted
     * @param hi the index of the end of the range to be sorted
     * @param cmp the comparison function to use
     */
    public def sortAtoms[T](a:Rail[T], lo:Long, hi:Long, cmp:(T,T)=>Int) {
        if (hi <= lo) return;
        var l:Long = lo - 1;
        var h:Long = hi;
        while (true) {
            while (cmp(a(++l),a(hi))<0);
            while (cmp(a(hi),a(--h))<0 && h>lo);
            if (l >= h) break;
            exch(a, l, h);
            swapAtoms(l, h);
        }
        exch(a, l, hi);
        swapAtoms(l, hi);
        sortAtoms[T](a, lo, l-1, cmp);
        sortAtoms[T](a, l+1, hi, cmp);
    }

    /**
     * Swap data between positions i and j in all of the particle data
     * arrays.
     */
    private @Inline final def swapAtoms(i:Long, j:Long):void {
        exch(globalIndex, i, j);
        exch(atomTypeIndex, i, j);
        exch(x, i, j);
        exch(fx, i, j);
    }

    private @Inline final def exch[T](a:Rail[T], i:Long, j:Long):void {
        val temp = a(i);
        a(i) = a(j);
        a(j) = temp;
    }

    private @Inline final def exch[T](a:ArrayList[T], i:Long, j:Long):void {
        val temp = a(i);
        a(i) = a(j);
        a(j) = temp;
    }
}

