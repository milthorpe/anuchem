/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2012.
 */
package au.edu.anu.mm;

import x10.compiler.Inline;
import x10.util.StringBuilder;

/**
 * This class calculates associated Legendre polynomials.
 * @author milthorpe
 */
public class AssociatedLegendrePolynomial {
    /** The terms P_{lm} (with m >= 0) in this polynomial */
    protected val terms:Rail[Double];

    /** The degree of the polynomial. */
    public val p:Long;

    public def this(p:Long) {
        this.terms = new Rail[Double]((p+1)*(p+1));
        this.p = p;
    }

    public @Inline final operator this(l:Long, m:Long) = terms(l*(l+1)+m);

    public @Inline final operator this(l:Long, m:Long)=(v:Double):Double {
        terms(l*(l+1)+m) = v;
        return v;
    }

    /**
     * Calculate associated Legendre polynomials P_{lk}(x) up to l=p (m.ge.0)
     * @see Dachsel (1996).
     *   "Fast and accurate determination of the Wigner rotation matrices in the fast multipole method".
     *   J. Chem. Phys. 124 (14) 144115. 14 April 2006.
     *   info:doi/10.1063/1.2194548
     */
    public static def getPlk(theta:Double, p:Long) : AssociatedLegendrePolynomial {
        val cosTheta = Math.cos(theta);
        val sinTheta = Math.sin(theta);

        val P = new AssociatedLegendrePolynomial(p);
        P(0,0) = 1.0;
        P(1,1) = -sinTheta;
        P(1,0) = cosTheta;
		var fact:Double = 1.0;
        for (l in 2..p) {
			fact += 2.0;
            P(l,l) = fact * sinTheta * -P(l-1,l-1);
            P(l,l-1) = fact * cosTheta * P(l-1,l-1);
            for (k in 0..(l-2)) {
            	P(l,k) = (fact * cosTheta * P(l-1,k) - (l+k-1) * P(l-2,k)) / (l-k);
            }
        }
        return P;
    }

    /**
     * Calculate associated Legendre polynomials P_{lk}(x) up to l=p (m.ge.0)
     * @see Dachsel (1996).
     *   "Fast and accurate determination of the Wigner rotation matrices in the fast multipole method".
     *   J. Chem. Phys. 124 (14) 144115. 14 April 2006.
     *   info:doi/10.1063/1.2194548
     */
    public def getPlk(theta:Double) {
        val cosTheta = Math.cos(theta);
        val sinTheta = Math.sin(theta);

        this(0,0) = 1.0;
        this(1,1) = -sinTheta;
        this(1,0) = cosTheta;
		var fact:Double = 1.0;
        for (l in 2..p) {
			fact += 2.0;
            this(l,l) = fact * sinTheta * -this(l-1,l-1);
            this(l,l-1) = fact * cosTheta * this(l-1,l-1);
            for (k in 0..(l-2)) {
            	this(l,k) = (fact * cosTheta * this(l-1,k) - (l+k-1) * this(l-2,k)) / (l-k);
            }
        }
    }

    /*
     * Calculate associated Legendre polynomials P_{lm}(x) up to l=p (m.ge.0)
     * @see White and Head-Gordon (1994).
     *      "Derivation and efficient implementation of the fast multipole method".
     *      J. Chem. Phys. 101 (8) 15 October 1994.
     */
    public static def getPlm(x:Double, p:Long) : AssociatedLegendrePolynomial {
        if (Math.abs(x) > 1.0) {
            throw new IllegalArgumentException("abs(x) > 1: Associated Legendre functions are only defined on [-1, 1].");
        }

        val Plm = new AssociatedLegendrePolynomial(p);
		//val Plm = DistArray.make[double](Region.makeLowerTriangular(p+1));
		Plm(0,0) = 1.0;
		val somx2 : Double = Math.sqrt((1.0 - x) * (1.0 + x));
		var fact : Double = 1.0;
		for (var i:Long = 1; i<=p; i++) {
			Plm(i,i) = -Plm(i-1,i-1) * fact * somx2;
			Plm(i,i-1) = x * fact * Plm(i-1,i-1);
			fact += 2.0;
		}
		for (var m:Long = 0; m<=p-2; m++) {
			for (var l:Long = m+2; l<=p; l++) {
				Plm(l,m) = (x * (2.0 * l - 1.0) 
							* Plm(l-1,m) - (l + m - 1.0) 
							* Plm(l-2,m)) / (l - m);
			}
		}
		return Plm;
	}

    /**
     * @return a string representation of this polynomial.  does not include terms for m<0
     */
    public def toString() : String {
        val s = new StringBuilder();
        for (i in 0..p) {
            for (j in 0..i) {
		        s.add("" + this(i,j) + " ");
            }
            s.add("\n");
	    }
        return s.toString();
    }
}
