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
    val placeId:Int;

    /** The number of atoms in all boxes below this box. */
    private var numAtoms:Int;

    /**
     * Creates a new GhostOctant for an octant at the given place.
     */
    public def this(id:OctantId, placeId:Int) {
        super(id, 0);
        this.placeId = placeId;
    }

    public def compareTo(b:Octant):Int = id.compareTo(b.id);

    public def countOctants():Int = 0;

    public def ghostOctants():Int = 1;

    public def numAtoms() = numAtoms;

    /** 
     * Go to home place of this octant and return multipole expansion (once computed).
     */
    protected def upward(localData:PlaceLocalHandle[FmmLocalData], size:Double, dMax:UByte):Pair[Int,MultipoleExpansion] {
        val mortonId = id.getMortonId();
        val result = at(Place(placeId)) getRemoteMultipole(localData, mortonId);
        numAtoms = result.first;
        //Console.OUT.println("at " + here + " GhostOctant.upward for " + id + " held at " + placeId + " numAtoms = " + numAtoms);
        return result;
    }

    private def getRemoteMultipole(localData:PlaceLocalHandle[FmmLocalData], mortonId:UInt):Pair[Int,MultipoleExpansion] {
        val octant = localData().getOctant(mortonId);
        if (octant != null) {
            //Console.OUT.println("at " + here + " waiting on multipole " + mortonId);
            when(octant.multipoleReady) {
                //Console.OUT.println("at " + here + " progressed on multipole " + mortonId + " numAtoms = " + octant.numAtoms);
                return Pair[Int,MultipoleExpansion](octant.numAtoms(), octant.multipoleExp);
            }
        } else {
            return Pair[Int,MultipoleExpansion](0, null);
        }
    }

    protected def downward(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, dMax:UByte):Double {
        //Console.OUT.println("at " + here + " GhostOctant.downward for " + id + " held at " + placeId); 
        if (numAtoms > 0) {
            return at(Place(placeId)) downwardRemote(localData, size, parentLocalExpansion, dMax);
        }
        return 0.0;
    }


    private def downwardRemote(localData:PlaceLocalHandle[FmmLocalData], size:Double, parentLocalExpansion:LocalExpansion, dMax:UByte):Double {
        val octant = localData().getOctant(id);
        if (octant != null) {
            return octant.downward(localData, size, parentLocalExpansion, dMax);
        } else {
            return 0.0;
        }
    }
     

    public def addToCombinedVSet(combinedVSet:HashSet[UInt], ws:Int) {
        super.addToCombinedVSet(combinedVSet, ws);
    }

    public def toString(): String {
        return "GhostOctant " + id;
    }
}

