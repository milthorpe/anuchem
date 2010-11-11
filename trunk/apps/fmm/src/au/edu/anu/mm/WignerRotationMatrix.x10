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

/**
 * This class calculates the Wigner rotation matrix D^l_{km}.
 * @see Dachsel (1996).
 *   "Fast and accurate determination of the Wigner rotation matrices in the fast multipole method".
 *   J. Chem. Phys. 124 (14) 144115. 14 April 2006.
 *   info:doi/10.1063/1.2194548
 * @author milthorpe
 */
public class WignerRotationMatrix {

	/**
	 * Calculate Wigner rotation matrix D_{km}(theta, l) for k = -l..l
	 */
	public static def getDkm(theta: double, l : int) : Array[Double](2) {
        if (Math.abs(theta) > 2.0 * Math.PI) {
            throw new IllegalArgumentException("abs(x) > 2*PI: Wigner rotation matrix is only defined on [0..2*PI].");
        }

        val Plm = AssociatedLegendrePolynomial.getPlm(Math.cos(theta), l);
        Console.OUT.println("Plm:");
		for ([i] in 0..l) {
            for ([j] in 0..i) {
			    Console.OUT.print("" + Plm(i,j) + " ");
            }
            Console.OUT.println();
		}

        val D = new Array[Double]((-l..l) * (-l..l));

        if (theta == 0) {
            // Eq. 30
            for ([k] in -l..l) {
                D(k,k) = 1.0;
            }
        } else if (theta == Math.PI) {
            // Eq. 30
            for ([k] in -l..l) {
                D(k,-k) = Math.pow(-1, l+k);
            }
        }

        var thetaPrime : Double = theta;
        // cosTheta < 0; need to calculate Dlm using addition theorems - Eq. 30
        if (theta > Math.PI / 2) {
            if (theta < Math.PI) {
                thetaPrime = Math.PI - theta;
            } else if (theta > Math.PI) {
                thetaPrime = theta - Math.PI;
            }
            // theta == PI already handled above
        }

        val cosTheta = Math.cos(thetaPrime);
        val sinTheta = Math.sin(thetaPrime);

        // gl0 starting point, Eq. 29
        var gk0 : Double = 1.0;
        for ([k] in 1..l) {
            gk0 = Math.sqrt((2.0*k-1) / (2.0*k)) * gk0;
        }

        // starting point for recursion, Eq. 28
        val gl0 = gk0;
        Console.OUT.println("gl0 = " + gl0);
        D(0,l) = Math.pow(-1.0, l) * gl0 * Math.pow(sinTheta, l);
        var glm : Double = gl0;
        var sign : Double = Math.pow(-1, l);
        for ([m] in 1..l) {
            glm = Math.sqrt((l-m+1) as Double / (l+m)) * glm; // Eq. 29
            sign *= -1.0;
            D(m,l) = sign * glm * Math.pow(1.0 + cosTheta, m) * Math.pow(sinTheta, l-m);
        }

        // Eq. 26
        for (var k:Int=l; k>-1; k--) {
            D(l,k-1) = (l+k) / Math.sqrt(l*(l+1.0) - k*(k-1.0)) * sinTheta / (1.0 + cosTheta) * D(l,k);
        }

        // Eq. 25
        for (var m:Int=l-1; m>=0; m--) {
            for (var k:Int=l; k>-l; k--) {
                D(m,k-1) = Math.sqrt(((l*(l+1) - m*(m+1)) as Double) / (l*(l+1) - k*(k-1))) * D(m+1,k)
                         + (m+k) / Math.sqrt((l*(l+1) - k*(k-1)) as Double) * sinTheta / (1.0 + cosTheta) * D(m,k);
            }
        }

        // Eq. 27
        for ([m] in -l..0) {
            sign = Math.pow(-1, m-l);
            for ([k] in -l..l) {
                D(m,k) = sign * D(-m,-k);
                sign *= -1;
            }
        }

		return D;
	}

}
