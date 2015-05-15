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
import x10.util.Pair;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * Holds all particle data for a given place in a molecular mechanics
 * simulation.  Atom features (globalIndex, atomType, residueNumber, position,
 * velocity, force) are held in separate Rails, indexed by a 'local atom index'
 * at this place. The globalIndex for an atom is unique across all places.
 */
public class ParticleData {
    public var description:String;

    public var atomTypes:Rail[AtomType];

    /** Global index of each atom */
    public val globalIndex = new ArrayList[Long]();
    /** Atom type index for each atom */
    public val atomTypeIndex = new ArrayList[Int]();
    /** Residue number (e.g. index of molecule) for each atom */
    public val residueNumber = new ArrayList[Int]();
    /** Atom centers in nm */
    public val x = new ArrayList[Point3d]();
    /** Atom velocities */
    public val dx = new ArrayList[Vector3d]();
    /** Instantaneous force on atoms */
    public val fx = new ArrayList[Vector3d]();

    public final def numAtoms() = globalIndex.size();

    /**
     * Exclusion list for non-bonded interactions. Each excluded pair is
     * entered twice as (a,b) and (b,a); the list is sorted by first atom
     * index to support searching.
     */
    protected val exclusions = new ArrayList[Pair[Long,Long]]();

    public compareExclusions = (a:Pair[Long,Long],b:Pair[Long,Long])=> {
        val f = a.first.compareTo(b.first);
        (f != 0n) ? f : a.second.compareTo(b.second)
    };

    public final @Inline def isExcluded(atomIndex1:Long, atomIndex2:Long) {
        val index = exclusions.binarySearch(Pair[Long,Long](atomIndex1,atomIndex2), compareExclusions);
        val excluded = (index >= 0 && exclusions(index).first == atomIndex1 && exclusions(index).second == atomIndex2);
        return excluded;
    }

    public def addAtom(index:Long, atomType:Int, center:Point3d, velocity:Vector3d) {
        globalIndex.add(index);
        atomTypeIndex.add(atomType);
        residueNumber.add(0n);
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
        residueNumber.removeAt(index);
        x.removeAt(index);
        dx.removeAt(index);
        fx.removeAt(index);
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
        if (lo < hi) {
            val p = partition(a, lo, hi, cmp);
            sortAtoms[T](a, lo, p-1, cmp);
            sortAtoms[T](a, p+1, hi, cmp);
        }
    }

    private @Inline def partition[T](a:Rail[T], lo:Long, hi:Long, cmp:(T,T)=>Int) {
        val pivotIndex = medianIndex(a, lo, hi, cmp);
        val pivotValue = a(pivotIndex);
        swapAtoms(a, pivotIndex, hi);

        var p:Long = lo;
        for (var i:Long = lo; i < hi; i++) {
            if (cmp(a(i), pivotValue) <= 0)
                swapAtoms(a, i, p++);
        }
        swapAtoms(a, p, hi);

        return p;
    }

    /**
     * Given indices lo and hi, returns the index of the median element of
     * a[lo, hi, (lo+hi)/2].
     */
    private @Inline def medianIndex[T](a:Rail[T], lo:Long, hi:Long, cmp:(T,T)=>Int) {
        val mid = (lo + hi) / 2;
        if (cmp(a(lo), a(mid)) < 0) {
            if (cmp(a(mid), a(hi)) < 0) {
                return mid;
            } else if (cmp(a(lo), a(hi)) < 0) {
                return hi;
            } else {
                return lo;
            }
        } else {
            if (cmp(a(lo), a(hi)) < 0) {
                return lo;
            } else if (cmp(a(mid), a(hi)) < 0) {
                return hi;
            } else {
                return mid;
            }
        }
    }

    /**
     * Swap data between positions i and j in all of the particle data
     * arrays.
     */
    private @Inline def swapAtoms[T](a:Rail[T], i:Long, j:Long):void {
        exch(a, i, j);
        exch(globalIndex, i, j);
        exch(atomTypeIndex, i, j);
        exch(residueNumber, i, j);
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

