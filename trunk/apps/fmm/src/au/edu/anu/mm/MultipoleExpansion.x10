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
 * This class calculates multipole expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public class MultipoleExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }


    public def this(p : Int, data : ValRail[Complex]) {
        super(p, data);
    }

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static safe def getOlm(q : Double, v : Tuple3d, p : Int) : MultipoleExpansion! {
        val exp = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Double = 1.0;
        for ((l) in 0..p) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac / ilm * (q * rfac * pplm(l,0)); 
            for ((m) in 1..l) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac / ilm * (q * rfac * pplm(l,m));
            }
            for ((m) in -l..-1) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate() * (2*((-m+1)%2)-1);
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /**
     * Calculate the chargeless multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static safe def getOlm(v : Tuple3d, p : Int) : MultipoleExpansion! {
        val exp = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Double = 1.0;
        for ((l) in 0..p) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac / ilm * (rfac * pplm(l,0)); 
            for ((m) in 1..l) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac / ilm * (rfac * pplm(l,m));
            }
            for ((m) in -l..-1) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate() * (2*((-m+1)%2)-1);
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /** 
     * Translate a multipole expansion centred around the origin along a vector -b,
     * and adds to this expansion centred at -b.
     * This corresponds to "Operator A", Equations 10-11 in White & Head-Gordon.
     * Note: this defines A^lm_jk(b) = O_l-j,m-k(b); however this is only defined
     * where abs(m-k) <= l-j, therefore we restrict m to [j-l+k..-j+l+k]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     */
    public def translateAndAddMultipole(shift : MultipoleExpansion!,
                                         source : MultipoleExpansion) {
        val p : Int = terms.region.max(0);
        val localSource = MultipoleExpansion.getLocalCopy(p, source);
        // TODO this atomic should be around inner loop update.
        // however as it's "stop the world" it's more efficient to do it out here
        atomic {
            for ((j,k) in terms.region) {
                val O_jk = localSource.terms(j,k);
                for ((l) in j..p) {
                    for ((m) in -l..l) {
                        if (Math.abs(m-k) <= (l-j)) {
                            val A_lmjk = shift.terms(l-j, m-k);
                            this.terms(l,m) = this.terms(l,m) + A_lmjk * O_jk;
                        }
                    }
                }
            }
        }
    }

    /**
     * For a periodic FMM, gets the macroscopic parent multipole
     * for the 3x3x3 box that is the parent of the box for which
     * this is the multipole expansion.
     * Uses the self similarity relation
     * M^(j+1)_(lm) = M^j_(lm) * 3^l
     * @see Kudin & Scuseria (1998) eq. 2.5
     */ 
    public def getMacroscopicParent() : MultipoleExpansion! {
        val p : Int = terms.region.max(0);
        val parentExpansion = new MultipoleExpansion(p);
        for ((l,m) in terms.region) {
            parentExpansion.terms(l,m) = terms(l,m) * Math.pow(3, l);
        }
        return parentExpansion;
    }

    /**
     * TODO this should not be necessary once XTENLANG-787 is resolved
     * @return a local copy of a multipole expansion from another place
     */
    static def getLocalCopy(p : Int, source : MultipoleExpansion) : MultipoleExpansion! {
        val data = at (source) {getData(p, source)};
        return new MultipoleExpansion(p, data);
    }
}


