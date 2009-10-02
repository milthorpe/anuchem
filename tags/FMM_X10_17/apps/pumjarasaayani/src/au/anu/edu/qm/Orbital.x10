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
    val exps:ArrayList[Double];
    val coeff:ArrayList[Double];
    val type:String;

    public def this(type:String) { 
       exps      = new ArrayList[Double]();
       coeff     = new ArrayList[Double]();
       this.type = type;
    } 

    public def addExponent(ex:Double) {
       exps.add(ex);
    }

    public def addCoefficient(co:Double) {
       coeff.add(co);
    }

    public def add(ex:Double, co:Double) {
       exps.add(ex);
       coeff.add(co);
    }

    public def getType() = this.type;
    public def getExponents() = this.exps;
    public def getCoefficients() = this.coeff;
}

