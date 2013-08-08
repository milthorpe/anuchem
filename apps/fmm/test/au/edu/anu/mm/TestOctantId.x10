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

import x10x.vector.Point3d;

/**
 * Test OctantId class
 * @author milthorpe
 */
public class TestOctantId {
    public static def main(args: Rail[String]) {
        val dMax = 5UY;
        val a = new OctantId(12UY,  2UY, 60UY, dMax-1UY);
        val b = new OctantId(11UY, 47UY, 21UY, dMax);
        val c = new OctantId(21UY, 30UY,  2UY, dMax);
        val d = new OctantId(12UY,  2UY, 60UY, dMax);
        Console.OUT.println(a.compareTo(b));
        Console.OUT.println(b.compareTo(c));
        Console.OUT.println(a.compareTo(c));
        Console.OUT.println(a.compareTo(a));
        Console.OUT.println(a.compareTo(d));

        Console.OUT.println("a = " + a + " parent = " + a.getParentId() + " offset = " + a.getParentId().getChildIndex(dMax, a));
        Console.OUT.println("b = " + b + " parent = " + b.getParentId() + " offset = " + b.getParentId().getChildIndex(dMax, b));
        Console.OUT.println("c = " + c + " parent = " + c.getParentId() + " offset = " + c.getParentId().getChildIndex(dMax, c));

        val x = new Point3d(0.0, 0.0, 0.0);
    }
}
