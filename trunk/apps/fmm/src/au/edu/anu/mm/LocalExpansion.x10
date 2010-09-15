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

import x10x.vector.Tuple3d;
import x10x.polar.Polar3d;

/**
 * This class calculates local Taylor-type expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public class LocalExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }

    /**
     * Calculate the local Taylor-type expansion M_{lm} (with m >= 0) for a single point v.
     */
    public static safe def getMlm(v : Tuple3d, p : int) : LocalExpansion {
        val exp = new LocalExpansion(p);
        val terms = exp.terms;
        val v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 

        val rfac0 : Double = 1.0 / v_pole.r;
        val phifac0 = Complex(Math.cos(v_pole.phi), Math.sin(v_pole.phi));
        var rfac : Double = rfac0;
        var il : Double = 1.0;
        for ([l] in 0..p) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            terms(l,0) = phifac * (rfac * pplm(l,0) * ilm);
            for ([m] in 1..l) {
                ilm = ilm / (l+1-m);
                phifac = phifac * phifac0;
                terms(l,m) = phifac * (rfac * pplm(l,m) * ilm);
            }
            for ([m] in -l..-1) {
                terms(l,m) = terms(l,-m).conjugate() * ((2*((-m+1)%2)-1 as Double));
            }
            rfac = rfac * rfac0;
        }

        return exp;
    }

    /** 
     * Translate a local expansion centred around the origin along a vector b,
     * and adds to this expansion centred at b, where shift is the expansion
     * of the vector b.
     * This corresponds to "Operator C", Equations 18-21 in White & Head-Gordon.
     * Note: this defines C^lm_jk(b) = O_j-l,k-m(b); however this is only defined
     * where abs(k-m) <= j-l, therefore we restrict k to [l-j+m..-l+j+m]
     * @param shift the multipole expansion of the translation
     * @param source the source local expansion, centred at the origin
     */
    public safe def translateAndAddLocal(shift : MultipoleExpansion,
                                         source : LocalExpansion) {
        val p = terms.region.max(0);

        // BEGIN HAND-INLINED ITERATOR
        // should be just:  for ([l,m] in terms) {
        val it_p = (terms.region as ExpansionRegion).p;
        var it_l:int = 0;
        var it_m:int = 0;
	while ((it_l <= it_p && it_m <= it_l)) {
            val l = it_l;
	    val m = it_m;
	    if (it_m<it_l) it_m++;
	    else {
                it_l++;
                it_m = -it_l;
            }
        // END HAND-INLINED ITERATOR
            for ([j] in l..p) {
                for ([k] in (l-j+m)..(-l+j+m)) {
                    val C_lmjk = shift.terms(j-l, k-m);
                    val O_jk = source.terms(j,k);
                    this.terms(l,m) = this.terms(l,m) + C_lmjk * O_jk;
                }
            }
        }
    }

   /** 
     * Transform a multipole expansion centred around the origin into a
     * Taylor expansion centred about b, and adds to this expansion.
     * This corresponds to "Operator B", Equations 13-15 in White & Head-Gordon.
     * This operator is inexact due to truncation of the series at <em>p</em> poles.
     * Note: this defines B^lm_jk(b) = M_j+l,k+m(b), therefore restrict l to [0..p-j]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     */
    public safe def transformAndAddToLocal(transform : LocalExpansion,
                                         source : MultipoleExpansion) {
        val p : Int = terms.region.max(0);

        // BEGIN HAND-INLINED ITERATOR
        // should be just:  for ([j,k] in terms) {
        val it_p = (terms.region as ExpansionRegion).p;
        var it_l:int = 0;
        var it_m:int = 0;
	while ((it_l <= it_p && it_m <= it_l)) {
            val j = it_l;
	    val k = it_m;
	    if (it_m<it_l) it_m++;
	    else {
                it_l++;
                it_m = -it_l;
            }
        // END HAND-INLINED ITERATOR
            val O_jk = source.terms(j,k);
            for ([l] in 0..p-j) {
                for ([m] in -l..l) {
                    if (Math.abs(k+m) <= (j+l)) {
                        val B_lmjk = transform.terms(j+l, k+m);
                        //Console.OUT.println("source.terms.dist(" + j + "," + k + ") = " + source.terms.dist(j,k));
                        this.terms(l,m) = this.terms(l,m) + B_lmjk * O_jk;
                    }
                }
            }
        }
    }

    /**
     * Transforms this local expansion about the origin to the potential
     * acting on <code>q</code> at point <code>v</code>.
     * @param q the charge at point v
     * @param v the location of charge q
     */
    public safe def getPotential(q : Double,
                                v : Tuple3d) : Double {
        val numTerms = terms.region.max(0);
        val transform = MultipoleExpansion.getOlm(q, v, numTerms);
        var potential : Double = 0.0;
        // TODO use lift/reduction?
        // BEGIN HAND-INLINED ITERATOR
        // should be just:  for ([i,j] in terms.region) {
        val it_p = (terms.region as ExpansionRegion).p;
        var it_l:int = 0;
        var it_m:int = 0;
	while ((it_l <= it_p && it_m <= it_l)) {
            val i = it_l;
	    val j = it_m;
	    if (it_m<it_l) it_m++;
	    else {
                it_l++;
                it_m = -it_l;
            }
        // END HAND-INLINED ITERATOR
            potential += (terms(i,j) * transform.terms(i,j)).re;
        }
        return potential;
    }

    /**
     * For a periodic FMM, gets the macroscopic parent local expansion
     * for the 3x3x3 box that is the parent of the box for which
     * this is the local expansion.
     * Uses the self similarity relation
     * L^(j+1)_(lm) = L^j_(lm) / 3^l+1
     * @see Kudin & Scuseria (1998) eq. 2.5
     */ 
    public def getMacroscopicParent() : LocalExpansion {
        val p : Int = terms.region.max(0);
        val parentExpansion = new LocalExpansion(p);
        for ([l,m] in terms.region) {
            parentExpansion.terms(l,m) = terms(l,m) / Math.pow(3.0, l+1);
        }
        return parentExpansion;
    }
}


