/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
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
		for ((i) : Point in [0..3]) {
            for ((j) : Point in [0..i]) {
			    Console.OUT.print(Plm(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestAssociatedLegendrePolynomial().execute();
    }

}
