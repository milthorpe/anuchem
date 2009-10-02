/**
 * ContractedGaussian.x10
 *
 * Represents a contracted gaussian function (made up of PrimitiveGaussian)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class ContractedGaussian { 
    val centeredAtom:Atom;
    val power:Power;
    var normalization:Double;
    val primitives:ArrayList[PrimitiveGaussian];
 
    val exponents:ArrayList[Double];
    val coefficients:ArrayList[Double];
    val primNorms:ArrayList[Double];

    public def this(atm:Atom, pwr:Power) { 
        this.centeredAtom = atm;
        this.power = pwr;
        
        normalization = 1.0; 
        primitives = new ArrayList[PrimitiveGaussian]();

        exponents = new ArrayList[Double]();
        coefficients = new ArrayList[Double]();
        primNorms = new ArrayList[Double]();
    } 

    public def getOrigin() = centeredAtom;
    public def getCenteredAtom() = centeredAtom;
    public def getPower() = power;
    public def getNormalization() = normalization;
    public def getPrimitives() = primitives;
    public def getExponents() = exponents;
    public def getCoefficients() = coefficients;
    public def getPrimNorms() = primNorms;
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    var index:Int;
    public def getIndex() = index;
    
    public def setIndex(idx:Int) {
       index = idx;
    }

    public def addPrimitive(exp:Double, coeff:Double) {
        primitives.add(new PrimitiveGaussian(centeredAtom, power, exp, coeff));

        exponents.add(exp);
        coefficients.add(coeff);
    }

    public def overlap(cg:ContractedGaussian) {
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

    public def kinetic(cg:ContractedGaussian) {
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

    public def nuclear(cg:ContractedGaussian, center:Atom) {
        var vij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();

        // TODO: x10 - parallel        
        for(i=0; i<primitives.size(); i++) {
            val iPG:PrimitiveGaussian = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                val jPG:PrimitiveGaussian = cgPrimitives.get(j);                
                
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

