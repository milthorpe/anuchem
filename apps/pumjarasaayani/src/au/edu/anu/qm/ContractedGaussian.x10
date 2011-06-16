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

import x10.compiler.NonEscaping;
import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * Represents a contracted Gaussian function (made up of PrimitiveGaussian)
 *
 * @author: V.Ganesh, milthorpe
 */
public struct ContractedGaussian { 
    public val centre : Point3d;
    public val power : Power;
    public val primitives : Rail[PrimitiveGaussian];
    public val normalization : Double;
    public val intIndex : Int;

    /**
     * Creates a new ContractedGaussian.
     * @param coeff the coefficients of the contracted Gaussian
     * @param exps the exponents of the contracted Gaussian.  Must be the same length as the coefficients array.
     */
    public def this(centre:Point3d, pwr:Power, primitives:Rail[PrimitiveGaussian], intIndex:Int, normalize:Boolean) { 
        this.centre = centre;
        this.power = pwr;
        this.primitives = primitives;
        this.intIndex = intIndex;
        if (normalize) {
            normalization = 1.0 / Math.sqrt(selfOverlap());
        } else {
            normalization = 1.0; 
        }
    } 

    public def getPrimitives() : Rail[PrimitiveGaussian] = primitives;
    public def getPrimitive(i:Int) : PrimitiveGaussian = primitives(i);
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    public def distanceFrom(cg:ContractedGaussian) : Double = centre.distance(cg.centre);
    public def distanceSquaredFrom(cg:ContractedGaussian) : Double = centre.distanceSquared(cg.centre);

    public def overlap(cg:ContractedGaussian) : Double {
        val cgPrimitives = cg.getPrimitives();
        var i:Int, j:Int;
        var sij:Double = 0.0;

        // TODO: x10 - parallel 
        for(i=0; i<primitives.size; i++) {
            val iPG = primitives(i);            
            for(j=0; j<cgPrimitives.size; j++) {
                val jPG = cgPrimitives(j);                
                
                sij += iPG.coefficient * jPG.coefficient * iPG.overlap(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * sij;
    }

    private @NonEscaping def selfOverlap() : Double {
        var sij:Double = 0.0;

        for(i in 0..(primitives.size-1)) {
            val iPG = primitives(i);            
            for(j in 0..(primitives.size-1)) {
                val jPG = primitives(j);                
                
                sij += iPG.coefficient * jPG.coefficient * iPG.overlap(jPG);
            }
        }
        
        return sij;
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
                
                tij += iPG.coefficient * jPG.coefficient * iPG.kinetic(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * tij;
    }

    public def nuclear(cg:ContractedGaussian, centre:Point3d) : Double {
        var vij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();

        // TODO: x10 - parallel        
        for(i=0; i<primitives.size; i++) {
            val iPG = primitives(i);            
            for(j=0; j<cgPrimitives.size; j++) {
                val jPG = cgPrimitives(j);                
                
                vij += iPG.coefficient * jPG.coefficient 
                       * iPG.nuclear(jPG, centre);                                
            } // end for
        } // end for
        
        return normalization * cg.normalization * vij;
    }
}

