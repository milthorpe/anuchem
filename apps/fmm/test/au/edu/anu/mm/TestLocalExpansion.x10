/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2014.
 */
package au.edu.anu.mm;

import x10x.vector.Vector3d;

/**
 * Test local Taylor-type expansions.
 * @author milthorpe
 */
class TestLocalExpansion extends MathTest {
    public def testAll(): boolean {
        val p = 3; // multipole level
        val x = Vector3d(1.0, 2.0, -1.0);

        val Mlm = LocalExpansion.getMlm(x, p);
        Console.OUT.println("local expansion:\n" + Mlm.toString());

        val target = new LocalExpansion(p);
        val y = Vector3d(2.0, -3.0, 1.0);
        target.translateAndAddLocal(y, Mlm);
        Console.OUT.println("translated local:\n" + target.toString());

        val roundtrip = new LocalExpansion(p);
        val z = Vector3d(-2.0, 3.0, -1.0);
        roundtrip.translateAndAddLocal(z, target);
        Console.OUT.println("translated local - roundtrip:\n" + roundtrip.toString());
		for (i in 0..p) {
            for (j in -i..i) {
                assert(nearlyEqual(roundtrip(i,j), Mlm(i,j), 1.0e-6, 1.0e-12));
		    }
        }

        return true;
    }

    public static def main(Rail[String]) {
        new TestLocalExpansion().testAll();
    }

}
