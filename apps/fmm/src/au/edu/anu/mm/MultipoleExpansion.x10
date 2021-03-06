/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010-2013.
 */
package au.edu.anu.mm;

import x10.regionarray.Array;

import x10x.vector.Vector3d;
import x10x.polar.Polar3d;

/**
 * This class calculates multipole expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public class MultipoleExpansion extends Expansion {

    public def this(p:Long) {
        super(p);
    }

    /**
     * A constructor which makes a copy of an existing expansion
     */
    public def this(source : MultipoleExpansion) {
	    super(source);
    }

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(q:Double, v:Vector3d, p:Long) : MultipoleExpansion {
        val exp = new MultipoleExpansion(p);
        val v_pole = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlk(v_pole.theta, p); 

        exp(0,0) = Complex(q, 0.0);

        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));
        var rfac : Double = v_pole.r;
        var il : Double = 1.0;
        for (l in 1..p) {
            il = il * l;
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp(l,0) = phifac / ilm * (q * rfac * pplm(l,0)); 
            for (m in 1..l) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
		        val O_lm = phifac / ilm * (q * rfac * pplm(l,m));
                exp(l,m) = O_lm;
    	    }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v
     * and add to this expansion.
     * @param q the charge
     * @param v the location of the charge relative to the expansion centre
     * @param plm a scratch space in which to store the assoc. Legendre polynomial for the centre
     */
    public def addOlm(q:Double, v:Vector3d, pplm:AssociatedLegendrePolynomial) {
        val v_pole = Polar3d.getPolar3d(v);
        pplm.getPlk(v_pole.theta); 

        this(0,0) += Complex(q, 0.0);

        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));
        var rfac : Double = v_pole.r;
        var il : Double = 1.0;
        for (l in 1..p) {
            il = il * l;
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            this(l,0) += phifac / ilm * (q * rfac * pplm(l,0)); 
            for (m in 1..l) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
		        val O_lm = phifac / ilm * (q * rfac * pplm(l,m));
                this(l,m) += O_lm;
    	    }
            rfac = rfac * v_pole.r;
        }
    }

    /**
     * Calculate the chargeless multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(v:Vector3d, p:Long) : MultipoleExpansion {
        val exp = new MultipoleExpansion(p);
        val v_pole = Polar3d.getPolar3d(v);
        val pplm = AssociatedLegendrePolynomial.getPlk(v_pole.theta, p); 

        exp(0,0) = Complex(1.0, 0.0);

        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));
        var rfac : Double = v_pole.r;
        var il : Double = 1.0;
        for (l in 1..p) {
            il = il * l;
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp(l,0) = phifac / ilm * (rfac * pplm(l,0)); 
            for (m in 1..l) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
        		val O_lm = phifac / ilm * (rfac * pplm(l,m));
                exp(l,m) = O_lm;
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /**
     * This is Operator A implementing rotations so that the actual translation occurs parallel with the z-axis
     * @param scratch, a MultipoleExpansion in which to perform temporary calculations
     * @param temp, a Complex Rail in which to store temporary results
     * @param v is the vector through which the source should be translated
     * @param complexK is the pre calculated values of exp(i*-k*phi)
     * @param source is the multipole to add
     * @param wigner is the pre calculated Wigner matrices for the rotation angle theta, indexed first by forwards (0) and backwards (1)
     * @see Dachsel 2006, eqn 9
     */
    public def translateAndAddMultipole(scratch:MultipoleExpansion, temp:Rail[Complex], v:Vector3d, complexK:Rail[Rail[Complex]], source:MultipoleExpansion, wigner:Rail[Rail[Array[Double](2){rect}]]) { 
	    val b = v.length();
	    val invB = 1 / b;

	    source.rotate(temp, complexK(0), wigner(0), scratch);

	    for (m in 0..p) {
		    for (l in m..p) temp(l) = scratch(l, m);

            var b_lm_pow:Double = 1.0;
		    for (l in m..p) {
			    var O_lm : Complex = Complex.ZERO;
			    var F_lm : Double = b_lm_pow / Factorial.getFactorial(l - m); // Factorial are already computed
			    for (j in m..l) {
				    O_lm = O_lm + temp(j) * F_lm; // explicitly this would be * Math.pow(translationPolar.r, l - j) / fact(l - j);
				    F_lm = F_lm * (l - j) * invB;
			    }
			    scratch(l, m) = O_lm;
                b_lm_pow = b_lm_pow * b;
		    }
	    }

	    scratch.backRotateAndAdd(complexK(1), wigner(1), this); 
    }

    /**
     * This is Operator A implementing rotations used if values are not already precalculated (i.e. in one-off tests)
     * Wigner matrices and values of exp(i*-k*phi) are calculated inside the method
     * @param v is the vector through which the source should be translated
     * @param source is the multipole to add
     */
    public def translateAndAddMultipole(v : Vector3d, source : MultipoleExpansion) {
        val scratch = new MultipoleExpansion(p);
        val temp = new Rail[Complex](p+1);
    	val polar = Polar3d.getPolar3d(v);
    	translateAndAddMultipole(scratch, temp, v, genComplexK(polar.phi, p), source, WignerRotationMatrix.getACollection(polar.theta, p) );
    }

    /**
     * Rotation of this expansion where Wigner matrices and exp(i*-k*phi) are not precalculated
     * Different method call for rotate which does the precalculations for the user
     * @param theta, rotation angle
     * @param phi, rotation angle
     * @return a new expansion which has been rotated
     */
    public def rotate(theta : Double, phi : Double) {
    	val target = new MultipoleExpansion(p);
        val temp = new Rail[Complex](p+1);
    	this.rotate(temp, genComplexK(phi, p)(1) , WignerRotationMatrix.getACollection(theta, p)(0), target);
    	return target;
    }

    /**
     * For a periodic FMM, gets the macroscopic parent multipole
     * for the 3x3x3 box that is the parent of the box for which
     * this is the multipole expansion.
     * Uses the self similarity relation
     * M^(j+1)_(lm) = M^j_(lm) * 3^l
     * @see Kudin & Scuseria (1998) eq. 2.5
     */ 
    public def getMacroscopicParent() : MultipoleExpansion {
        val parentExpansion = new MultipoleExpansion(p);
        for (l in 0..p) {
            for (m in 0..l) {
                parentExpansion(l,m) = this(l,m) * Math.pow(3.0, l);
            }
        }
        return parentExpansion;
    }
}
