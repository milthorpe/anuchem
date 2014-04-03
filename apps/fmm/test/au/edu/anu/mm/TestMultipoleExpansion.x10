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
 * Test multipole expansions
 * @author milthorpe
 */
class TestMultipoleExpansion extends MathTest {
    public def testAll(): boolean {
        val p = 3; // multipole level
        val x = Vector3d(1.0, 2.0, -1.0);
        val Olm = MultipoleExpansion.getOlm(1.5, x, p);

        Console.OUT.println("multipole expansion:\n" + Olm.toString());

        val target = new MultipoleExpansion(p);
        val translation = MultipoleExpansion.getOlm(Vector3d(2.0, -3.0, 1.0), p);
        target.translateAndAddMultipole(translation, Olm);
        Console.OUT.println("translated multipole:\n" + target.toString());

        val roundtrip = new MultipoleExpansion(p);
        val reverseTranslation = MultipoleExpansion.getOlm(Vector3d(-2.0, 3.0, -1.0), p);
        roundtrip.translateAndAddMultipole(reverseTranslation, target);
        Console.OUT.println("translated multipole - roundtrip:\n" + roundtrip.toString());
		for (i in 0..p) {
            for (j in -i..i) {
                assert(nearlyEqual(roundtrip(i,j), Olm(i,j), 1.0e-6, 1.0e-12)); 
		    }
        }

        val localExp = new LocalExpansion(p);
        val transformation = LocalExpansion.getMlm(Vector3d(2.0, -3.0, 1.0), p);
        localExp.transformAndAddToLocal(transformation, Olm);
        Console.OUT.println("transformed multipole:\n" + localExp.toString());

        val rotated = Olm.rotate(0.0/*Math.PI/4.0*/, 3.0*Math.PI/4.0);
        Console.OUT.println("rotated multipole:\n" + rotated.toString());
        val restored = rotated.rotate(0.0/*-Math.PI/4.0*/, -3.0*Math.PI/4.0);
        Console.OUT.println("restored multipole:\n" + restored.toString());

        return true;
    }

    public static def main(Rail[String]) {
        new TestMultipoleExpansion().testAll();
    }
}
