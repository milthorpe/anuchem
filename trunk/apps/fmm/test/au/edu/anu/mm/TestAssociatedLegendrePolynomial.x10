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

import harness.x10Test;
import au.edu.anu.mm.AssociatedLegendrePolynomial;

/**
 * Test calculation of associated Legendre polynomials.
 * @author milthorpe
 */
class TestAssociatedLegendrePolynomial extends x10Test {
    public def run(): boolean {
        val Plm = AssociatedLegendrePolynomial.getPlm(0.5, 3);
		for ((i) : Point in 0..3) {
            for ((j) : Point in 0..i) {
			    Console.OUT.print(Plm(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Array[String](1)) {
        new TestAssociatedLegendrePolynomial().execute();
    }

}
