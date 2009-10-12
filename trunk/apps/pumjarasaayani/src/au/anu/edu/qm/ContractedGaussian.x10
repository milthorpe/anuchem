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
    global var centeredAtom:Atom{self.at(this)};
    global var power:Power{self.at(this)};
    global var normalization:Double;
    global var primitives:ArrayList[PrimitiveGaussian{self.at(this)}]{self.at(this)};
 
    global var exponents:ArrayList[Double]{self.at(this)};
    global var coefficients:ArrayList[Double]{self.at(this)};
    global var primNorms:ArrayList[Double]{self.at(this)};

    public def this() { }

    public def make(atm:Atom{self.at(this)}, pwr:Power{self.at(this)}) { 
        this.centeredAtom = atm;
        this.power = pwr;
        
        normalization = 1.0; 
        primitives = new ArrayList[PrimitiveGaussian{self.at(this)}]();

        exponents = new ArrayList[Double]();
        coefficients = new ArrayList[Double]();
        primNorms = new ArrayList[Double]();
    } 

    public def getOrigin() : Atom{self.at(this)} = centeredAtom;
    public def getCenteredAtom() : Atom{self.at(this)} = centeredAtom;
    public def getPower() : Power{self.at(this)} = power;
    public def getNormalization() = normalization;
    public def getPrimitives() : ArrayList[PrimitiveGaussian{self.at(this)}]{self.at(this)} = primitives;
    public def getExponents() : ArrayList[Double]{self.at(this)} = exponents;
    public def getCoefficients() : ArrayList[Double]{self.at(this)} = coefficients;
    public def getPrimNorms() : ArrayList[Double]{self.at(this)} = primNorms;
    public def getTotalAngularMomentum() = power.getTotalAngularMomentum();
    public def getMaximumAngularMomentum() = power.getMaximumAngularMomentum();
    public def getMinimumAngularMomentum() = power.getMinimumAngularMomentum();

    var index:Int;
    public def getIndex() = index;
    
    public def setIndex(idx:Int) {
       index = idx;
    }

    public def addPrimitive(exp:Double, coeff:Double) {
        val pg:PrimitiveGaussian{self.at(this)} = new PrimitiveGaussian();
        pg.make(centeredAtom, power, exp, coeff);
        primitives.add(pg);

        exponents.add(exp);
        coefficients.add(coeff);
    }

    public def overlap(cg:ContractedGaussian{self.at(this)}) : Double {
        val cgPrimitives = cg.getPrimitives();
        var i:Int, j:Int;
        var sij:Double = 0.0;

        // TODO: x10 - parallel 
        for(i=0; i<primitives.size(); i++) {
            var iPG:PrimitiveGaussian{self.at(this)} = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                var jPG:PrimitiveGaussian{self.at(this)} = cgPrimitives.get(j);                
                
                sij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.overlap(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * sij;
    }

    public def kinetic(cg:ContractedGaussian{self.at(this)}) : Double {
        var tij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();
        
        // TODO: x10 - parallel 
        for(i=0; i<primitives.size(); i++) {
            var iPG:PrimitiveGaussian{self.at(this)} = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                var jPG:PrimitiveGaussian{self.at(this)} = cgPrimitives.get(j);                
                
                tij += iPG.getCoefficient() * jPG.getCoefficient() * iPG.kinetic(jPG);
            } // end for
        } // end for
        
        return normalization * cg.normalization * tij;
    }

    public def nuclear(cg:ContractedGaussian{self.at(this)}, center:Atom{self.at(this)}) : Double {
        var vij:Double = 0.0;
        var i:Int, j:Int;
        val cgPrimitives = cg.getPrimitives();

        // TODO: x10 - parallel        
        for(i=0; i<primitives.size(); i++) {
            val iPG:PrimitiveGaussian{self.at(this)} = primitives.get(i);            
            for(j=0; j<cgPrimitives.size(); j++) {
                val jPG:PrimitiveGaussian{self.at(this)} = cgPrimitives.get(j);                
                
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

