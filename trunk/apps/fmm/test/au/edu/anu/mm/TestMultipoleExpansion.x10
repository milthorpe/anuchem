/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Tuple3d;

/**
 * Test multipole expansions
 * @author milthorpe
 */
class TestMultipoleExpansion extends MathTest {
    public def run(): boolean {
        val p = 3; // multipole level

        val x = new Point3d(1.0, 2.0, -1.0);
        val Olm = MultipoleExpansion.getOlm(1.5, x as Tuple3d, p);
        Console.OUT.println("multipole expansion:\n" + Olm.toString());

        val target = new MultipoleExpansion(p);
        val translation = MultipoleExpansion.getOlm(new Point3d(2.0, -3.0, 1.0) as Tuple3d, p);
        target.translateAndAddMultipole(translation, Olm);
        Console.OUT.println("translated multipole:\n" + target.toString());

        val roundtrip = new MultipoleExpansion(p);
        val reverseTranslation = MultipoleExpansion.getOlm(new Point3d(-2.0, 3.0, -1.0) as Tuple3d, p);
        roundtrip.translateAndAddMultipole(reverseTranslation, target);
        Console.OUT.println("translated multipole - roundtrip:\n" + roundtrip.toString());
		for ([i] in 0..p) {
            for ([j] in -i..i) {
                chk(nearlyEqual(roundtrip.terms(i,j), Olm.terms(i,j), 1.0e-6, 1.0e-12)); 
		    }
        }

        val localExp = new LocalExpansion(p);
        val transformation = LocalExpansion.getMlm(new Point3d(2.0, -3.0, 1.0) as Tuple3d, p);
        localExp.transformAndAddToLocal(transformation, Olm);
        Console.OUT.println("transformed multipole:\n" + localExp.toString());

        // test copy
        val localCopy = MultipoleExpansion.getLocalCopy(p, roundtrip);
		for ((i): Point in 0..p) {
            for ((j): Point in -i..i) {
                chk(roundtrip.terms(i,j) == localCopy.terms(i,j));
		    }
        }        

        return true;
    }

    public static def main(Array[String](1)) {
        new TestMultipoleExpansion().execute();
    }

}
