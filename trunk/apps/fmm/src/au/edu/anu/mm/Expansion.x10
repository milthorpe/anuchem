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

import x10.util.StringBuilder;

/**
 * This is the superclass for multipole and local expansions, as used in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * These expansions have a peculiarly shaped region(abs(x0)<=x1 && 0<=x1<=p)
 * (X10 gives it as (x0+x1>=0 && x0-x1>=0 && x0<=3), which is constructed
 * here by subtracting two halfspaces from a rectangular region. 
 * @author milthorpe
 */
public class Expansion {
    /** The terms X_{lm} (with m >= 0) in this expansion */
    public val terms : Array[Complex](2);
    public static val ones = new Array[Complex](-100..100, Complex.ONE);

    public def this(p : Int) {
        //var expRegion : Region(2) = [0..p,-p..p];
        //expRegion = expRegion - Region.makeHalfspace([1,1],1);
        //expRegion = expRegion - Region.makeHalfspace([1,-1],1);
        val expRegion = new ExpansionRegion(p);
        this.terms = new Array[Complex](expRegion);
    }

    public def this(e : Expansion) { 
	    this.terms = new Array[Complex](e.terms);
    }

    public atomic def add(e : Expansion) {
	val p = terms.region.max(0);
	    for ([l] in -p..p) {
	        for ([m] in -l..l) { // this will now be inlined -- should be for ([l,m] in terms.region)
	            this.terms(l,m) = this.terms(l,m) + e.terms(l,m);
    	    }
        }
    }

    public def toString() : String {
        val p : Int = terms.region.max(0);
        val s = new StringBuilder();
        for ([i] in 0..p) {
            for ([j] in -i..i) {
		        s.add("" + terms(i,j) + " ");
            }
            s.add("\n");
	    }
        return s.toString();
    }

    /**
     * Code to calculate the exponentials ( exp(i*k*phi) ) required to do a rotation around z-axis given an angle, used in Fmm3d and also for test cases where once-off angles are needed
     * @param phi, k
     * @return array of Complex first indexed by forward (0), backward (1) then by k
     * @see Dachsel 2006, eqn 4 & 5
     */
    public static def genComplexK(phi : Double, p : int) { 
    	val complexK = new Array[Array[Complex](1)](0..1);
    	for ([r] in 0..1) { 
	    	complexK(r) = new Array[Complex](-p..p); 
    		for ([k] in -p..p) complexK(r)(k) = Math.exp(Complex.I * k * phi * ((r==0)?1:-1) );
	    }
    	return complexK;
    }
    public static def genComplexKZero(p : int) = ones; // quicker code for when angle is 0 

    /** 
     * Rotates this expansion (local and multipole are differentiated by different precalculated wigner matrices) in three dimensions.
     * @param wigner, precalculated wigner matrices, an array of WignerMatrices, indexed first by p from 0 to max terms in expansion
     * @param complexK, values of exp(i*k*phi)
     * @see Dachsel 2006 eqn 4 & 5
     */
    public def rotate( complexK : Array[Complex](1), wigner : Array[Array[Double](2){rect}](1) ) {
        val p : Int = terms.region.max(0);

    	val lkFac = new Array[Complex](-p..p);
    	var O_lm : Complex;
        for ([l] in 1..p) {
            val Dl = wigner(l); // avoids calculating matrices directly

	        for ([k] in -l..l) lkFac(k) = terms(l, k) * complexK(k);
           
	        var m_sign : int = 1;
            for ([m] in 0..l) {
	            O_lm = Complex.ZERO;
                for ([k] in -l..l) {
                    O_lm = O_lm + lkFac(k) * Dl(m, k); // Eq. 5
                }
                terms(l,m) = O_lm;

        	    if (m != 0) terms(l, -m) = O_lm.conjugate() * m_sign;
            	m_sign = -m_sign; // instead of doing the conjugate
            }
        }
    }

    /**
     * Special case of expansion rotation when theta = 0, to use when unrotating an expansion after translation (also no return value)
     * @param complexK, same meaning as above
     */
    public def phiRotate(complexK : Array[Complex](1) ) {
	val p = terms.region.max(0);
        for ([l] in -p..p) {
            for ([m] in -l..l) {
                terms(l, m) = complexK(m) * terms(l, m);
            }
        }
    }

}
