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

public struct OctantId(x:UByte, y:UByte, z:UByte, level:UByte) implements Comparable[OctantId] {
    static TOP_LEVEL = 2UY;

    public def this(x:UByte, y:UByte, z:UByte, level:UByte) {
        property(x,y,z,level);
    }

    public def compareTo(b:OctantId):Int {
        return this.getMortonId().compareTo(b.getMortonId());
    }

    public operator this < (x:OctantId) = (this.compareTo(x) < 0);
    public operator this <= (x:OctantId) = (this.compareTo(x) <= 0);
    public operator this > (x:OctantId) = (this.compareTo(x) > 0);
    public operator this >= (x:OctantId) = (this.compareTo(x) >= 0);

    public def next():OctantId {
        return OctantId.getFromMortonId(this.getMortonId()+(1U << 8));
    }

    public def getMortonId():UInt {
        val x = this.x as UInt;
        val y = this.y as UInt;
        val z = this.z as UInt;
        //Console.OUT.printf("x %8s y %8s z %8s level %8s\n", x.toBinaryString(), y.toBinaryString(), z.toBinaryString(), level.toBinaryString());
        var id:UInt = 0U;
        var bitmask:UInt = 1U;
        var shift:Int = 8;
        for (i in 0..7) {
            id |= (bitmask & z) << shift++;
            id |= (bitmask & y) << shift++;
            id |= (bitmask & x) << shift;
            bitmask = bitmask << 1;
        }
        id |= (level as UInt);
        //Console.OUT.printf("Morton id = %32s\n", id.toBinaryString());
        return id;
    }

    public static def getFromMortonId(mortonId:UInt):OctantId {
        //Console.OUT.printf("get from Morton id = %32s\n", mortonId.toBinaryString());
        val level = mortonId & UByte.MAX_VALUE as UInt;
        var x:UInt = 0U;
        var y:UInt = 0U;
        var z:UInt = 0U;
        var shift:Int = 24;
        var bitmask:UInt = 1U << 31;
        for (i in 0..7) {
            x |= (mortonId & bitmask) >> shift--; bitmask = bitmask >> 1;
            y |= (mortonId & bitmask) >> shift--; bitmask = bitmask >> 1;
            z |= (mortonId & bitmask) >> shift;   bitmask = bitmask >> 1;
        }
        //Console.OUT.printf("x %8s y %8s z %8s level %8s\n", (x as UByte).toBinaryString(), (y as UByte).toBinaryString(), (z as UByte).toBinaryString(), (level as UByte).toBinaryString());
        return OctantId(x as UByte, y as UByte, z as UByte, level as UByte);
    }
        

    public def getLeafMortonId():UInt {
        val x = this.x as UInt;
        val y = this.y as UInt;
        val z = this.z as UInt;
        //Console.OUT.printf("x %8s y %8s z %8s\n", x.toBinaryString(), y.toBinaryString(), z.toBinaryString());
        var id:UInt = 0U;
        var bitmask:UInt = 1U;
        var shift:Int = 0;
        for (i in 0..7) {
            id |= (bitmask & z) << shift++;
            id |= (bitmask & y) << shift++;
            id |= (bitmask & x) << shift;
            bitmask = bitmask << 1;
        }
        //Console.OUT.printf("leaf Morton id = %32s\n", id.toBinaryString());
        return id;
    }

    /** @return the octant ID of the parent of this octant */
    public def getParentId():OctantId {
        return new OctantId(x/2, y/2, z/2, level-1);
    }

    /** @return the anchor of this octant */
    public def getAnchor(dMax:UByte):OctantId {
        if (level == dMax) return this;
        val d = (Math.pow2(dMax) as UByte) / (Math.pow2(level) as UByte); 
        return new OctantId(x*d, y*d, z*d, dMax);
    }

    /** @return the index of the given child octant in this (shared) octant's child array */
    public def getChildIndex(dMax:UByte, childId:OctantId):Int {
        return (childId.x%2 << 2) | (childId.y%2 << 1) | (childId.z%2);
    }

    public def toString(): String {
        return "" + level + ":(" + x + "," + y + "," + z + ")";
    }
}
