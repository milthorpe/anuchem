/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */
package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.polar.Polar3d;
import x10x.vector.Vector3d;

/**
 * This class is a container for the Wigner rotation matrices 
 * and exponential prefactors required to perform translation
 * and transformation operations on multipole and local
 * expansions in the Fast Multipole Method.
 * @see Dachsel (1996).
 *   "Fast and accurate determination of the Wigner rotation matrices in the fast multipole method".
 *   J. Chem. Phys. 124 (14) 144115. 14 April 2006.
 *   info:doi/10.1063/1.2194548
 * @author milthorpe
 */
class FmmOperators {
    /** 
     * A cache of wigner matrices for translating between box centres of children to parent
     * Dimensions of the Array are:
     *   0, 1, 2: x, y, z translation values
     * which in turn contains an Array indexed by:
     *   0 for forward, 1 for backward (rotations)
     * Each element of this is an Array[Array[Double](2){rect}](1)  
     *   which is indexed by l-value
     * Each element of this is a Wigner rotation matrix d^l for a particular  
     *   theta (for the translation with vector <x,y,z>)
     */
    public val wignerA : Array[Rail[Rail[Array[Double](2){rect}]]](3){rect,zeroBased};
    public val wignerB : Array[Rail[Rail[Array[Double](2){rect}]]](3){rect}; // B is not zero-based
    public val wignerC : Array[Rail[Rail[Array[Double](2){rect}]]](3){rect,zeroBased};

    /**
     * A cache of exp(k * phi * i) values for all phi that could be needed in a rotation and -p < k < p
     * Dimensions of the Array are:
     *   0, 1, 2: x, y, z translation values
     * which in turn contains an Array indexed by:
     *   0 for +phi and 1 for -phi (for forward, back rotations)
     */
    public val complexK : Array[Rail[Array[Complex](1){rect,rail==false}]](3){rect};

    /**
     * Initialises the FMM operator arrays.
     * @param numTerms number of terms in multipole and local expansions
     * @param ws well-separated parameter
     */
    public def this(numTerms:Int, ws:Int) {
    	this.wignerA = precomputeWignerA(numTerms);
    	this.wignerB = precomputeWignerB(numTerms, ws);
    	this.wignerC = precomputeWignerC(numTerms);
    	this.complexK = precomputeComplex(numTerms, ws);
    }

/** 
     * Precomputes wigner rotation matrices premultiplied by appropriate factors for use 
     * in applying operator A. This is replicated in a unique dist across all places.
     */
    private def precomputeWignerA(numTerms : Int) {
        val wignerA = new Array[Rail[Rail[Array[Double](2){rect}]]]((0..1)*(0..1)*(0..1));
        for ([i,j,k] in wignerA) {
    		val theta = Polar3d.getPolar3d( Point3d(i*2-1,j*2-1,k*2-1) ).theta;
	        wignerA(i, j, k) = WignerRotationMatrix.getACollection(theta, numTerms);
        }
        return wignerA;
    }

    private def precomputeWignerC(numTerms : Int) {
        val wignerB = new Array[Rail[Rail[Array[Double](2){rect}]]]((0..1)*(0..1)*(0..1));
        for ([i,j,k] in wignerB) {
		    val theta = Polar3d.getPolar3d( Point3d(i*2-1,j*2-1,k*2-1) ).theta;
	        wignerB(i, j, k) = WignerRotationMatrix.getCCollection(theta, numTerms);
        }
        return wignerB;
    }

    private def precomputeWignerB(numTerms : Int, ws : Int) {
        val wignerC = new Array[Rail[Rail[Array[Double](2){rect}]]]((-(2*ws+1))..(2*ws+1) * (-(2*ws+1))..(2*ws+1) * (-(2*ws+1))..(2*ws+1));
        for ([i,j,k] in wignerC) {
            val theta = Polar3d.getPolar3d ( Point3d(i, j, k) ).theta;
	        wignerC(i, j, k) = WignerRotationMatrix.getBCollection(theta, numTerms);
        }
        return wignerC;
    }

    /**
     * Precomputes the values of exp(phi * k * i) needed for every possible translation operation
     * and replicates to all places using a unique dist
     */
    private def precomputeComplex(numTerms : Int, ws : Int) {
        val complexK = new Array[Rail[Array[Complex](1){rect,rail==false}]]((-(2*ws+1))..(2*ws+1) * (-(2*ws+1))..(2*ws+1) * (-(2*ws+1))..(2*ws+1));
        for ([i,j,k] in complexK) {
            val phi = Polar3d.getPolar3d ( Point3d(i, j, k) ).phi;
            complexK(i, j, k) = Expansion.genComplexK(phi, numTerms);
        }
        return complexK;
    }

}
