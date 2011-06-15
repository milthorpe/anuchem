/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * Represents a contracted gaussian function (made up of PrimitiveGaussian)
 *
 * @author: V.Ganesh, milthorpe
 */
public class ContractedGaussian { 
    val center : Point3d;
    val power : Power;
    val primitives : Rail[PrimitiveGaussian];

    var normalization : Double;

    /**
     * Creates a new ContractedGaussian.
     * @param coeff the coefficients of the contracted Gaussian
     * @param exps the exponents of the contracted Gaussian.  Must be the same length as the coefficients array.
     */
    public def this(center:Point3d, pwr:Power, coeff:ArrayList[Double], exps:ArrayList[Double]) { 
        this.center = center;
        this.power = pwr;
        normalization = 1.0; 

        primitives = new Array[PrimitiveGaussian](coeff.size());
        for(var i:Int=0; i<coeff.size(); i++) {
            primitives(i) = new PrimitiveGaussian(center, power, exps(i), coeff(i));
        }
    } 

    public def getOrigin() = center;
    public def getCenter() = center;
    public def getPower() = power;
    public def getNormalization() = normalization;
    public def setNormalization(n:Double) : void { normalization = n; }
    public def getPrimitives() : Rail[PrimitiveGaussian] = primitives;
    public def getPrimitive(i:Int) : PrimitiveGaussian = primitives(i);
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    var intIndex:Int;
    public def getIntIndex() = intIndex;

    public def setIntIndex(idx:Int) {
       intIndex = idx;
    }

    public def distanceFrom(cg:ContractedGaussian) : Double = center.distance(cg.center);
    public def distanceSquaredFrom(cg:ContractedGaussian) : Double = center.distanceSquared(cg.center);

    public def overlap(cg:ContractedGaussian) : Double {
        val cgPrimitives = cg.getPrimitives();
        var i:Int, j:Int;
        var sij:Double = 0.0;

        // TODO: x10 - parallel 
        for(i=0; i<primitives.size; i++) {
            val iPG = primitives(i);            
            for(j=0; j<cgPrimitives.size; j++) {
                val jPG = cgPrimitives(j);                
                
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
        for(i=0; i<primitives.size; i++) {
            val iPG = primitives(i);            
            for(j=0; j<cgPrimitives.size; j++) {
                val jPG = cgPrimitives(j);                
                
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
        for(i=0; i<primitives.size; i++) {
            val iPG = primitives(i);            
            for(j=0; j<cgPrimitives.size; j++) {
                val jPG = cgPrimitives(j);                
                
                vij += iPG.getCoefficient() * jPG.getCoefficient() 
                       * iPG.nuclear(jPG, center);                                
            } // end for
        } // end for
        
        return normalization * cg.normalization * vij;
    }

    public def normalize() {        
        for(var i:Int=0; i<primitives.size; i++) {
           primitives(i).normalize();
        }

        normalization = 1.0 / Math.sqrt(this.overlap(this));
    }
}

