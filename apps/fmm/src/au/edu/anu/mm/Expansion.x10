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

import x10.array.TriangularRegion;

import x10.compiler.Inline;
import x10.util.StringBuilder;

/**
 * This is the superclass for multipole and local expansions, as used in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * These expansions have a peculiarly shaped region(abs(x1)<=x0 && 0<=x0<=p)
 * however only terms with x1>=0 are stored, as the terms with x1<0 are just
 * complex conjugates of the terms with x1>0.
 * @author milthorpe
 */
public class Expansion {
    /** The terms X_{lm} (with m >= 0) in this expansion */
    protected val terms : Rail[Complex];

    /** The number of terms in the expansion. */
    public val p : Int;

    public def this(p : Int) {
        this.terms = new Array[Complex]((p+1)*(p+1));
        this.p = p;
    }

    public def this(e : Expansion) { 
	    this.terms = new Array[Complex](e.terms);
        this.p = e.p;
    }

    public @Inline final operator this(l:Int, m:Int) = terms(l*(l+1)+m);

    public @Inline final operator this(l:Int, m:Int)=(v:Complex):Complex {
        terms(l*(l+1)+m) = v;
        return v;
    }

    public def clear() { terms.clear(); }

    /**
     * Add each term of e to this expansion. 
     * This operation is atomic and therefore thread-safe.
     */
    public atomic def add(e : Expansion) {
        unsafeAdd(e);
    }

    /**
     * Add each term of e to this expansion. 
     * This operation is not atomic, therefore not thread-safe.
     */
    @Inline def unsafeAdd(e : Expansion) {
        for (i in terms) {
            terms(i) += e.terms(i);
        }
    }

    /**
     * @return a string representation of this expansion.  does not include terms for m<0
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

    /**
     * Code to calculate the exponentials ( exp(i*k*phi) ) required to do a rotation around z-axis given an angle, used in Fmm3d and also for test cases where once-off angles are needed
     * @param phi, k
     * @return array of Complex first indexed by forward (0), backward (1) then by k
     * @see Dachsel 2006, eqn 4 & 5
     */
    public static def genComplexK(phi:Double, p:Int):Rail[Rail[Complex]] { 
        val complexK = new Array[Rail[Complex]](2);
    	for (r in 0..1) {
            complexK(r) = new Array[Complex](p+1);
            for (k in 0..p) complexK(r)(k) = Math.exp(Complex.I * k * phi * ((r==0)?1:-1) );
	    }
    	return complexK;
    }

    /** 
     * Rotates this expansion (local and multipole are differentiated by different precalculated wigner matrices) in three dimensions.
     * Performs rotation around z-axis first, THEN rotation around x-axis (for "forwards" rotation)
     * @param temp an Rail[Complex] of size at least p+1 in which to perform temporary calculations
     * @param wigner, precalculated wigner matrices, an array of WignerMatrices, indexed first by p from 0 to max terms in expansion
     * @param complexK, values of exp(i*k*phi)
     * @see Dachsel 2006 eqn 4 & 5
     */
    public def rotate(temp:Rail[Complex], complexK:Rail[Complex], wigner:Rail[Array[Double](2){rect}]) {
        for (l in 1..p) {
            val Dl = wigner(l);

	        for (k in 0..l) temp(k) = this(l, k) * complexK(k);
           
            for (m in 0..l) {
	            var O_lm:Complex = temp(0) * Dl(m, 0);
                var m_sign:Double = -1.0;
                for (k in 1..l) {
                    val temp_k = temp(k);
                    O_lm += temp_k * Dl(m, k) + m_sign * temp_k.conjugate() * Dl(m, -k); // Eq. 5, for k and -k
                    m_sign = -m_sign;
                }
                this(l,m) = O_lm;
            }
        }
    }

    /** 
     * Rotates this expansion (local and multipole are differentiated by different precalculated wigner matrices) in three dimensions.
     * Performs rotation around x-axis first, THEN rotation around z-axis (for "backwards" rotation)
     * @param temp an Rail[Complex] of size at least p+1 in which to perform temporary calculations
     * @param wigner, precalculated wigner matrices, an array of WignerMatrices, indexed first by p from 0 to max terms in expansion
     * @param complexK, values of exp(i*k*phi)
     * @see Dachsel 2006 eqn 4 & 5
     */
    public def backRotate(temp:Rail[Complex], complexK:Rail[Complex], wigner:Rail[Array[Double](2){rect}]) {
        for (l in 1..p) {
            val Dl = wigner(l);

	        for (k in 0..l) temp(k) = this(l, k);
           
            for (m in 0..l) {
	            var O_lm:Complex = temp(0) * Dl(m, 0);
                var m_sign:Double = -1.0;
                for (k in 1..l) {
                    val temp_k = temp(k);
                    O_lm += temp_k * Dl(m, k) + m_sign * temp_k.conjugate() * Dl(m, -k); // Eq. 5, for k and -k
                    m_sign = -m_sign;
                }
                O_lm = O_lm * complexK(m);
                this(l,m) = O_lm;
            }
        }
    }

    /** 
     * Rotates this expansion (local and multipole are differentiated by different precalculated wigner matrices) in three dimensions and adds to the given target expansion
     * Performs rotation around x-axis first, THEN rotation around z-axis (for "backwards" rotation)
     * @param wigner, precalculated wigner matrices, an array of WignerMatrices, indexed first by p from 0 to max terms in expansion
     * @param complexK, values of exp(i*k*phi)
     * @see Dachsel 2006 eqn 4 & 5
     */
    public def backRotateAndAdd(complexK:Rail[Complex], wigner:Rail[Array[Double](2){rect}], target:Expansion) {
        target(0,0) += this(0,0);
        for (l in 1..p) {
            val Dl = wigner(l);
            val t_l0 = this(l,0);

            for (m in 0..l) {
	            var O_lm:Complex = t_l0 * Dl(m, 0);
                var m_sign:Double = -1.0;
                for (k in 1..l) {
                    val t_lk = this(l, k);
                    O_lm += t_lk * Dl(m, k) + m_sign * t_lk.conjugate() * Dl(m, -k); // Eq. 5, for k and -k
                    m_sign = -m_sign;
                }
                O_lm = O_lm * complexK(m);
                target(l,m) += O_lm;
            }
        }
    }

}
