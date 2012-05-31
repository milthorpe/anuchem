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

public struct OctantId(x:UShort, y:UShort, z:UShort, level:UShort) implements Comparable[OctantId] {
    public def this(x:UShort, y:UShort, z:UShort, level:UShort) {
        property(x,y,z,level);
    }

    public def compareTo(b:OctantId):Int {
        return this.getMortonId().compareTo(b.getMortonId());
    }

    public def getMortonId():ULong {
        val x = this.x as ULong;
        val y = this.y as ULong;
        val z = this.z as ULong;
        //Console.OUT.printf("x %8s y %8s z %8s level %8s\n", x.toBinaryString(), y.toBinaryString(), z.toBinaryString(), level.toBinaryString());
        var id:ULong = 0L;
        var bitmask:ULong = 1L;
        var shift:Int = 8;
        for (i in 0..7) {
            id |= (bitmask & z) << shift++;
            id |= (bitmask & y) << shift++;
            id |= (bitmask & x) << shift++;
            bitmask = bitmask << 1;
        }
        id |= (level as ULong);
        //Console.OUT.printf("Morton id = %64s\n", id.toBinaryString());
        return id;
    }

    /** @return the anchor of the parent of this octant */
    public def getParentId(dMax:UShort):OctantId {
        val d = (Math.pow2(dMax) as UShort) / (Math.pow2(level-1) as UShort); 
        return new OctantId(x-x%d, y-y%d, z-z%d, level-1);
    }

    /** @return the index of the given child octant in this (shared) octant's child array */
    public def getChildIndex(dMax:UShort, childId:OctantId):Int {
        val d = (Math.pow2(dMax) as UShort) / (Math.pow2(childId.level) as UShort); 
        return (((childId.x/d)%2) << 2) | (((childId.y/d)%2) << 1) | ((childId.z/d)%2);
    }

    public def toString(): String {
        return "" + level + ":(" + x + "," + y + "," + z + ")";
    }
}
