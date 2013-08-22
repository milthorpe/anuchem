/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */
package x10x.vector;

import harness.x10Test;

/**
 * Tests XLA Vector class
 * @author milthorpe
 */
public class TestVector extends x10Test {
    public def run(): boolean {
        val vector = new Vector(3);
        for(i in 0..(vector.vec.size-1)) {
            vector.vec(i) = i as Double;
        }
        val magnitude : Double = vector.magnitude();
        chk(magnitude == Math.sqrt(5.0));
        return true;
    }

    public static def main(Rail[String]) {
        new TestVector().execute();
    }
}
