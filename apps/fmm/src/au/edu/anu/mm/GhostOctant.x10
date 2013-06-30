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

/**
 * This class represents a ghost octant, i.e. an octant for which data are held
 * at another place, in the 3D division of space for the Fast Multipole Method.
 * @author milthorpe
 */
public class GhostOctant extends Octant implements Comparable[Octant] {
    /** The host place of the target octant for which this GhostOctant is a proxy. */
    val placeId:Long;

    /** The target octant. */
    private var target:GlobalRef[Octant];

    /** The number of atoms in all boxes below this box. */
    private var numAtoms:Long;

    /**
     * Creates a new GhostOctant for an octant at the given place.
     */
    public def this(id:OctantId, placeId:Long) {
        super(id);
        this.placeId = placeId;
    }

    public def compareTo(b:Octant):Int = id.compareTo(b.id);

    public def countOctants():Int = 0;

    public def ghostOctants():Int = 1;

    public def numAtoms() = numAtoms;

    static class GhostUpward(
        numAtoms:Long,
        multipoleExp:MultipoleExpansion,
        target:GlobalRef[Octant]
    ) {

        public def this(numAtoms:Long, multipoleExp:MultipoleExpansion, target:GlobalRef[Octant]) {
            property(numAtoms, multipoleExp, target);
        }
    }

    /** 
     * Go to home place of this octant and return multipole expansion (once computed).
     */
    protected def upward():Pair[Long,MultipoleExpansion] {
        val mortonId = id.getMortonId();
        val result = at(Place(placeId)) GhostOctant.upwardRemote(mortonId);
        if (result != null) {
            numAtoms = result.numAtoms;
            target = result.target;
            //Console.OUT.println("at " + here + " GhostOctant.upward for " + id + " held at " + placeId + " numAtoms = " + numAtoms);
            return Pair[Long,MultipoleExpansion](result.numAtoms, result.multipoleExp);
        } else {
            return Pair[Long,MultipoleExpansion](0L, null);
        }
    }

    private static def upwardRemote(mortonId:UInt):GhostUpward {
        val octant = FastMultipoleMethod.localData.getOctant(mortonId);
        if (octant != null) {
            //Console.OUT.println("at " + here + " waiting on multipole " + mortonId);
            when(octant.multipoleReady);
            //Console.OUT.println("at " + here + " progressed on multipole " + mortonId + " numAtoms = " + octant.numAtoms);
            return new GhostUpward(octant.numAtoms(), octant.multipoleExp, new GlobalRef[Octant](octant));
        } else {
            return null;
        }
    }

    protected def downward(parentLocalExpansion:LocalExpansion):Double {
        //Console.OUT.println("at " + here + " GhostOctant.downward for " + id + " held at " + placeId); 
        if (numAtoms > 0) {
            val target = this.target;
            return at (target.home) target().downward(parentLocalExpansion);
        }
        return 0.0;
    }

    public def toString(): String {
        return "GhostOctant " + id;
    }
}

