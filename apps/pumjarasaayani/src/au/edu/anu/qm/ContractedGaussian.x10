/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.compiler.NonEscaping;
import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * Represents a contracted Gaussian function (made up of PrimitiveGaussian)
 *
 * @author: V.Ganesh, milthorpe
 */
public struct ContractedGaussian { 
    public val origin : Point3d;
    public val power : Power;
    public val normalization : Double;
    public val exponents:Rail[Double];
    public val coefficients:Rail[Double];
    public val intIndex : Int;

    /**
     * Creates a new ContractedGaussian.
     * @param coeff the coefficients of the contracted Gaussian
     * @param exps the exponents of the contracted Gaussian.  Must be the same length as the coefficients array.
     */
    public def this(origin:Point3d, pwr:Power, exponents:Rail[Double], coefficients:Rail[Double], intIndex:Int, normalize:Boolean) { 
        this.origin = origin;
        this.power = pwr;
        this.exponents = exponents;
        this.coefficients = coefficients;
        this.intIndex = intIndex;
        if (normalize) {
            normalization = 1.0 / Math.sqrt(selfOverlap());
        } else {
            normalization = 1.0; 
        }
    } 

    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMaximumDegreeOfContraction() = coefficients.size as Int;

    public def distanceFrom(cg:ContractedGaussian) : Double = origin.distance(cg.origin);
    public def distanceSquaredFrom(cg:ContractedGaussian) : Double = origin.distanceSquared(cg.origin);

    public def overlap(cg:ContractedGaussian):Double {
        val cgExponents = cg.exponents;
        val cgCoefficients = cg.coefficients;
        val cgOrigin = cg.origin;
        val cgPower = cg.power;

        // TODO: x10 - parallel
        var sij:Double = 0.0;
        for(var i:Long=0; i<exponents.size; i++) {
            val expI = exponents(i);
            val coeffI = coefficients(i);
            val normI = PrimitiveGaussian.getNormalization(power, expI);
            for(var j:Long=0; j<cgExponents.size; j++) {
                val expJ = cgExponents(j);
                val coeffJ = cgCoefficients(j);
                val normJ = PrimitiveGaussian.getNormalization(cgPower, expJ);
                sij += coeffI * coeffJ * normI * normJ * PrimitiveGaussian.overlap(expI, origin, power, expJ, cgOrigin, cgPower);
            } // end for
        } // end for
        
        return normalization * cg.normalization * sij;
    }

    private @NonEscaping def selfOverlap() : Double {
        var sij:Double = 0.0;

        for(i in 0..(exponents.size-1)) {
            val expI = exponents(i);
            val coeffI = coefficients(i);
            val normI = PrimitiveGaussian.getNormalization(power, expI);
            for(j in 0..(exponents.size-1)) {
                val expJ = exponents(j);
                val coeffJ = coefficients(j);
                val normJ = PrimitiveGaussian.getNormalization(power, expJ);
                sij += coeffI * coeffJ * normI * normJ * PrimitiveGaussian.overlap(expI, origin, power, expJ, origin, power);
            }
        }
        
        return sij;
    }

    public def kinetic(cg:ContractedGaussian) : Double {
        val cgExps = cg.exponents;
        val cgCoeffs = cg.coefficients;
        val cgOrigin = cg.origin;
        val cgPower = cg.power;

        // TODO: x10 - parallel 
        var tij:Double = 0.0;
        for(var i:Long=0; i<exponents.size; i++) {
            val expI = exponents(i);
            val coeffI = coefficients(i);
            val normI = PrimitiveGaussian.getNormalization(power, expI);
            for(var j:Long=0; j<cgExps.size; j++) {
                val expJ = cgExps(j);
                val coeffJ = cgCoeffs(j);
                val normJ = PrimitiveGaussian.getNormalization(cgPower, expJ);
                tij += coeffI * coeffJ * PrimitiveGaussian.kinetic(expI, origin, power, normI, expJ, cgOrigin, cgPower, normJ);
            } // end for
        } // end for
        
        return normalization * cg.normalization * tij;
    }

    public def nuclear(cg:ContractedGaussian, centre:Point3d) : Double {
        val cgExps = cg.exponents;
        val cgCoeffs = cg.coefficients;
        val cgOrigin = cg.origin;
        val cgPower = cg.power;

        // TODO: x10 - parallel        
        var vij:Double = 0.0;
        for(var i:Long=0; i<exponents.size; i++) {
            val iPG = new PrimitiveGaussian(origin, power, exponents(i), coefficients(i), true);
            for(var j:Long=0; j<cgExps.size; j++) {
                val jPG = new PrimitiveGaussian(cgOrigin, cgPower, cgExps(j), cgCoeffs(j), true);
                
                vij += iPG.coefficient * jPG.coefficient 
                       * iPG.nuclear(jPG, centre);                                
            } // end for
        } // end for
        
        return normalization * cg.normalization * vij;
    }
}

