/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * ContractedGaussian.x10
 *
 * Represents a contracted gaussian function (made up of PrimitiveGaussian)
 *
 * @author: V.Ganesh
 */
public class ContractedGaussian { 
    val center : Point3d;
    val power : Power;
    val primitives : ArrayList[PrimitiveGaussian];
    val exponents:ArrayList[Double];
    val coefficients:ArrayList[Double];
    val primNorms:ArrayList[Double];

    val maxam:Int, minam:Int, totam:Int;

    var normalization : Double;

    public def this(center:Point3d, pwr:Power) { 
        this.center = center;
        this.power = pwr;
        normalization = 1.0; 
        primitives = new ArrayList[PrimitiveGaussian]();

        exponents = new ArrayList[Double]();
        coefficients = new ArrayList[Double]();
        primNorms = new ArrayList[Double]();

        maxam = power.getMaximumAngularMomentum();
        minam = power.getMinimumAngularMomentum();
        totam = power.getTotalAngularMomentum();
    } 

    public def getOrigin() = center;
    public def getCenter() = center;
    public def getPower() = power;
    public def getNormalization() = normalization;
    public def setNormalization(n:Double) : Void { normalization = n; }
    public def getPrimitives() : ArrayList[PrimitiveGaussian] = primitives;
    public def getPrimitive(i:Int) : PrimitiveGaussian = primitives.get(i);
    public def getExponents() : ArrayList[Double] = exponents;
    public def getCoefficients() : ArrayList[Double] = coefficients;
    public def getPrimNorms() : ArrayList[Double] = primNorms;
    public def getPrimNorm(i:Int) = primNorms(i);
    public def getTotalAngularMomentum() = totam;
    public def getMaximumAngularMomentum() = maxam;
    public def getMinimumAngularMomentum() = minam;

    var index:Int;
    public def getIndex() = index;
    
    public def setIndex(idx:Int) {
       index = idx;
    }

    var intIndex:Int;
    public def getIntIndex() = intIndex;

    public def setIntIndex(idx:Int) {
       intIndex = idx;
    }

    public def distanceFrom(cg:ContractedGaussian) : Double = center.distance(cg.center);
    public def distanceSquaredFrom(cg:ContractedGaussian) : Double = center.distanceSquared(cg.center);

    public def addPrimitive(exp:Double, coeff:Double) {
        val pg = new PrimitiveGaussian(center, power, exp, coeff);
        primitives.add(pg);

        exponents.add(exp);
        coefficients.add(coeff);
    }

    public def overlap(cg:ContractedGaussian) : Double {
        val cgPrimitives = cg.getPrimitives();
        var i:Int, j:Int;
        var sij:Double = 0.0;

        // TODO: x10 - parallel 
        for(i=0; i<primitives.size(); i++) {
            var iPG:PrimitiveGaussian = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                var jPG:PrimitiveGaussian = cgPrimitives.get(j);                
                
                sij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.overlap(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * sij;
    }

    public def kinetic(cg:ContractedGaussian) : Double {
        var tij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();
        
        // TODO: x10 - parallel 
        for(i=0; i<primitives.size(); i++) {
            var iPG:PrimitiveGaussian = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                var jPG:PrimitiveGaussian = cgPrimitives.get(j);                
                
                tij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.kinetic(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * tij;
    }

    public def nuclear(cg:ContractedGaussian, center : Point3d) : Double {
        var vij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();

        // TODO: x10 - parallel        
        for(i=0; i<primitives.size(); i++) {
            val iPG = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                val jPG = cgPrimitives.get(j);                
                
                vij += iPG.getCoefficient() * jPG.getCoefficient() 
                       * iPG.nuclear(jPG, center);                                
            } // end for
        } // end for
        
        return normalization * cg.normalization * vij;
    }

    public def normalize() {        
        for(var i:Int=0; i<primitives.size(); i++) {
           primitives.get(i).normalize();
           primNorms.add(primitives.get(i).getNormalization());
        }

        normalization = 1.0 / Math.sqrt(this.overlap(this));
    }
}

