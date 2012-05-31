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

/**
 * Test OctantId class
 * @author milthorpe
 */
class TestOctantId {
    public static def main(Rail[String]) {
        val dMax = 5US;
        val a = new OctantId(12,  2, 60, dMax-1);
        val b = new OctantId(11, 47, 21, dMax);
        val c = new OctantId(21, 30,  2, dMax);
        Console.OUT.println(a.compareTo(b));
        Console.OUT.println(b.compareTo(c));
        Console.OUT.println(a.compareTo(c));

        Console.OUT.println("a = " + a + " parent = " + a.getParentId(dMax) + " offset = " + a.getParentId(dMax).getChildIndex(dMax, a));
        Console.OUT.println("b = " + b + " parent = " + b.getParentId(dMax) + " offset = " + b.getParentId(dMax).getChildIndex(dMax, b));
        Console.OUT.println("c = " + c + " parent = " + c.getParentId(dMax) + " offset = " + c.getParentId(dMax).getChildIndex(dMax, c));
    }
}
