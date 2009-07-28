/**
 * Orbitals.x10
 *
 * Represents an Orbital (used for basis function formation)
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class Orbital { 
    var exps:ArrayList[Double];
    var coeff:ArrayList[Double];
    var type:String;

    public def this(type:String) { 
       exps      = new ArrayList[Double]();
       coeff     = new ArrayList[Double]();
       this.type = type;
    } 

    public def addExponent(ex:Double) : void {
       exps.add(ex);
    }

    public def addCoefficient(co:Double) : void {
       coeff.add(co);
    }

    public def add(ex:Double, co:Double) : void {
       exps.add(ex);
       coeff.add(co);
    }

    public def getType() : String = this.type;
    public def getExponents() : ArrayList[Double] = this.exps;
    public def getCoefficients() : ArrayList[Double] = this.coeff;
}

