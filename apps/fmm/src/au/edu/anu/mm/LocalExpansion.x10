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

import x10.regionarray.Array;

import au.edu.anu.chem.mm.MMAtom;
import x10x.vector.Vector3d;
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
     * A constructor which makes a copy of an existing expansion
     */
    public def this(source : LocalExpansion) {
    	super(source);
    }

    /**
     * Calculate the local Taylor-type expansion M_{lm} (with m >= 0) for a single point v.
     */
    public static def getMlm(v:Vector3d, p:Int) : LocalExpansion {
        val exp = new LocalExpansion(p);
        val v_pole = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlk(v_pole.theta, p);
        val rfac0 : Double = 1.0 / v_pole.r;

        exp(0,0) = Complex(rfac0, 0.0);

        val phifac0 = Complex(Math.cos(v_pole.phi), Math.sin(v_pole.phi));
        var rfac : Double = rfac0 * rfac0;
        var il : Double = 1.0;
        for (l in 1..p) {
            il = il * l;
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp(l,0) = phifac * (rfac * pplm(l,0) * ilm);
            for (m in 1..l) {
                ilm = ilm / (l+1-m);
                phifac = phifac * phifac0;
        		val M_lm = phifac * (rfac * pplm(l,m) * ilm);
                exp(l,m) = M_lm;
            }
            rfac = rfac * rfac0;
        }

        return exp;
    }

    /** 
     * More efficient version of Operator C (translation and addition) with rotations
     * @param scratch, a MultipoleExpansion in which to perform temporary calculations
     * @param temp, a Complex Rail in which to store temporary results
     * @param v a Tuple representing the shift of the expansion
     * @param wigner a collection of Wigner matrices precalculated to speed up the rotation
     * @param complexK is the pre calculated values of exp(i*k*phi)
     * @param source the source local expansion
     * @see Dachsel 2006, eqn 11
     */
    public def translateAndAddLocal(scratch:MultipoleExpansion, temp:Rail[Complex], v:Vector3d, complexK:Rail[Rail[Complex]], source:LocalExpansion, wigner:Rail[Rail[Array[Double](2){rect}]]) { 
	    val b = v.length();

        source.rotate(temp, complexK(1), wigner(0), scratch);

    	for (m in 0..p) {
    		for (l in m..p) temp(l) = scratch(l, m);

    		for (l in m..p) {
    			var M_lm : Complex = Complex.ZERO;
    			var F_lm : Double = 1.0;
    			for (j in l..p) {  
    				M_lm = M_lm + temp(j) * F_lm;
    				F_lm = F_lm * b / (j - l + 1);
    			}
    			scratch(l, m) = M_lm;
    		}
	   	}

	    scratch.backRotateAndAdd(complexK(0), wigner(1), this);
    }

    /**
     * Different method call for Operator C which pre calculates matrices and exp(k*i*phi) for use in tests etc 
     * @param v a Tuple representing the shift of the expansion
     * @param source the source local expansion
     */
    public def translateAndAddLocal(v : Vector3d, source : LocalExpansion) {
        val scratch = new MultipoleExpansion(p);
        val temp = new Rail[Complex](p+1);
	    val polar = Polar3d.getPolar3d(v);
	    translateAndAddLocal(scratch, temp, v, genComplexK(polar.phi, p), source, WignerRotationMatrix.getCCollection(polar.theta, p) );
    }

    /** 
     * More efficient version of Operator B with rotations
     * @param scratch, a MultipoleExpansion which can be used to do temporary calculations in
     * @param temp, a Complex Rail to do temporary calculations in
     * @param v, a Tuple representing the shift of the expansion
     * @param complexK, is the pre calculated values of exp(i*k*phi)
     * @param source, the expansion which should be shifted and added to this one
     * @param wigner, a collection of Wigner matrices precalculated to speed up the rotation
     * @see Dachsel 2006, eqn 10
     */
    public def transformAndAddToLocal(scratch:MultipoleExpansion, temp:Rail[Complex], v:Vector3d, complexK:Rail[Rail[Complex]], source:MultipoleExpansion, wigner:Rail[Rail[Array[Double](2){rect}]]) { 
        val inv_b = 1 / v.length();

        source.rotate(temp, complexK(0), wigner(0), scratch);

        var m_sign:Double = 1.0;
        var b_m_pow:Double = 1.0;
        for (m in 0..p) {
            for (l in m..p) {
                temp(l) = m_sign * scratch(l, m).conjugate();
            }

            var F_lm:Double = Factorial.getFactorial(m+m) * inv_b * b_m_pow * b_m_pow;
            var ml:Double = m+m+1;
            for (l in m..p) {
                var M_lm : Complex = Complex.ZERO;
                var F_jl : Double = F_lm;
                var jl:Double = ml;

                for (j in m..p) {
                    M_lm = M_lm + temp(j) * F_jl;
                    F_jl = jl * inv_b * F_jl;
                    jl += 1.0;
                }
                scratch(l, m) = M_lm;
                F_lm = ml * inv_b * F_lm;
                ml += 1;
            }
            m_sign = -m_sign;
            b_m_pow = b_m_pow * inv_b;
        }

        scratch.backRotateAndAdd(complexK(0), wigner(1), this);
    }

    /**
     * Different method call for Operator B which pre calculates matrices and exp(k*i*phi) for use in tests etc 
     * @param v, a Tuple representing the shift of the expansion
     * @param source, the source local expansion which will be added to this one
     */
    public def transformAndAddToLocal(v : Vector3d, source : MultipoleExpansion) {
    	val polar = Polar3d.getPolar3d(v);
        val scratch = new MultipoleExpansion(p);
        val temp = new Rail[Complex](p+1);
    	transformAndAddToLocal(scratch, temp, v, genComplexK(polar.phi, p), source, WignerRotationMatrix.getBCollection(polar.theta, p) );
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
    public def transformAndAddToLocal(transform : LocalExpansion,
                                         source : MultipoleExpansion) {
        // TODO should be just:  for ([j,k] in terms.region) {
        for (j in 0..p) {
            var k_sign:Double=1-(2*j%2);
            for (k in -j..j) {
                val O_jk = k < 0 ? (k_sign * source(j,-k).conjugate()) : source(j,k);
                for (l in 0..(p-j)) {
                    for (m in 0..l) {
                        val km = k+m;
                        if (Math.abs(km) <= (j+l)) {
                            val B_lmjk = km < 0 ? ((1-(2*km%2)) * transform(j+l, -km).conjugate()) : transform(j+l, km);
                            //Console.OUT.println("source.terms.dist(" + j + "," + k + ") = " + source.terms.dist(j,k));
                            this(l,m) = this(l,m) + B_lmjk * O_jk;
                        }
                    }
                }
                k_sign = -k_sign;
            }
        }
    }

    /**
     * Calculate the potential and forces for an atom due to this expansion,
     * and add the forces to the atom.
     * @param atom the atom for which to calculate the potential
     * @param v the vector from the box centre to the atom centre
     * @param plm a scratch space in which to store the assoc. Legendre polynomial for the atom centre
     * @return the potential due to distant particles
     * @see Kabadshow (2006). "The Fast Multipole Method - Alternative Gradient Algorithm and Parallelization". PhD thesis, Forschungszentrum Juelich.  Section 3.2.1
     */
    public def calculatePotentialAndForces(atom:MMAtom, v:Vector3d, pplm:AssociatedLegendrePolynomial):Double {
        val v_pole = Polar3d.getPolar3d(v);
        val q = atom.charge;
        pplm.getPlk(v_pole.theta);

        var dr:Double = 0.0;
        var dt:Double = 0.0;
        var dp:Double = 0.0;

        var potential:Double = q * this(0,0).re;

        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));
        var rfac : Double = v_pole.r;
        var rfacPrev : Double = 1.0;
        var il : Double = 1.0;
        for (l in 1..p) {
            val Ml0 = this(l,0);
            il = il * l;
            var F_lm:Complex = Complex(1.0/il, 0.0); // e^{-i m phi} / (l+|m|)!
            potential += (Ml0 * F_lm * (q * rfac * pplm(l,0))).re;
            dr += (Ml0 * F_lm * (q * l * rfacPrev * pplm(l,0))).re;
            dt += (Ml0 * F_lm * (q * rfacPrev * -pplm(l,1))).re;
            // phi terms cancel for m=0
            for (m in 1..l) {
                F_lm = F_lm * phifac0 / (l+m);
                val Olm = F_lm * (q * rfac * pplm(l,m));
                val Mlm = this(l,m);

                potential += 2.0*(Mlm.re * Olm.re) - 2.0*(Mlm.im * Olm.im);

                val r_lm = F_lm * (q * l * rfacPrev * pplm(l,m));
                dr += 2.0*(Mlm.re * r_lm.re) - 2.0*(Mlm.im * r_lm.im); // avoids conjugate for mirror terms m < 0
                val Plm1 = (m<l) ? pplm(l,m+1) : 0.0;
                val theta_lm = F_lm * 0.5 * q * rfacPrev * ((l-m+1)*(l+m) * pplm(l,m-1) - Plm1);
                dt += 2.0*(Mlm.re * theta_lm.re) - 2.0*(Mlm.im * theta_lm.im); // avoids conjugate for mirror terms m < 0
                val phi_lm = Complex.I * F_lm * -0.5 * q * rfacPrev * ((l-m+1)*(l-m+2) * pplm(l+1,m-1) + pplm(l+1,m+1));
                dp += 2.0*(Mlm.re * phi_lm.re) - 2.0*(Mlm.im * phi_lm.im); // avoids conjugate for mirror terms m < 0
    	    }
            rfacPrev = rfac;
            rfac = rfac * v_pole.r;
        }
        atom.force += v_pole.getGradientVector(dr, -dt, -dp);
        return potential;
    }

    /**
     * Different method call for rotate which does the precalculations for the user
     * @param theta, rotation angle
     * @param phi, rotation angle
     * @return a new expansion, which is the current expansion after it is rotated
     */
    public def rotate(theta : Double, phi : Double) {
        val target = new LocalExpansion(p);
        val temp = new Rail[Complex](p+1);
    	this.rotate(temp, genComplexK(phi, p)(0), WignerRotationMatrix.getCCollection(theta, p)(0), target);
    	return target;
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
        val parentExpansion = new LocalExpansion(p);
        for (l in 0..p) {
            for (m in 0..l) {
                parentExpansion(l,m) = this(l,m) / Math.pow(3.0, l+1);
            }
        }
        return parentExpansion;
    }
}


