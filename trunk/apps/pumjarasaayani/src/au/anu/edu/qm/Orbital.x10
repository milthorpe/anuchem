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

import x10.util.*;

/**
 * Orbitals.x10
 *
 * Represents an Orbital (used for basis function formation)
 *
 * @author: V.Ganesh
 */
public class Orbital { 
    global val exps:ArrayList[Double]{self.at(this)};
    global val coeff:ArrayList[Double]{self.at(this)};
    global val type:String{self.at(this)};
    global val angularMomentum:Int;

    public def this(type:String) { 
       exps      = new ArrayList[Double]();
       coeff     = new ArrayList[Double]();
       this.type = type as String{self.at(this)};

       if (type.equals("S")) angularMomentum = 0;
       else if (type.equals("P")) angularMomentum = 1;
       else if (type.equals("D")) angularMomentum = 2;
       else if (type.equals("F")) angularMomentum = 3;
       else angularMomentum = 0;
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

    public def getType() : String{self.at(this)} = this.type;
    public def getAngularMomentum() = this.angularMomentum;
    public def getExponents() : ArrayList[Double]{self.at(this)} = this.exps;
    public def getCoefficients() : ArrayList[Double]{self.at(this)} = this.coeff;
}

