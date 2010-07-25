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

import x10.array.TriangularRegion;

/**
 * This class calculates associated Legendre polynomials.
 * @author milthorpe
 */
public class AssociatedLegendrePolynomial {

	/**
	 * Calculate associated Legendre polynomials P_{lm}(x) up to l=p (m.ge.0)
	 */
	public static safe def getPlm(x: double, p : int) : Array[Double](2)! {
        if (Math.abs(x) > 1.0) {
            throw new IllegalArgumentException("abs(x) > 1: Associated Legendre functions are only defined on [-1, 1].");
        }

        val triRegion = new TriangularRegion(0,0,p+1,true);
        val Plm = new Array[Double](triRegion);
		//val Plm = DistArray.make[double](Region.makeLowerTriangular(p+1));
		Plm(0,0) = 1.0;
		val somx2 : Double = Math.sqrt((1.0 - x) * (1.0 + x));
		var fact : Double = 1.0;
		for (var i : Int = 1; i<=p; i++) {
			Plm(i,i) = -Plm(i-1,i-1) * fact * somx2;
			Plm(i,i-1) = x * fact * Plm(i-1,i-1);
			fact += 2.0;
		}
		for (var m : Int = 0; m<=p-2; m++) {
			for (var l : Int = m+2; l<=p; l++) {
				Plm(l,m) = (x * (2.0 * l - 1.0) 
							* Plm(l-1,m) - (l + m - 1.0) 
							* Plm(l-2,m)) / (l - m);
			}
		}
		return Plm;
	}

}
