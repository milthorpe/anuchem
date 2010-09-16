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
    val exps:ArrayList[Double];
    val coeff:ArrayList[Double];
    val type:String;
    val angularMomentum:Int;

    public def this(type:String) { 
       exps      = new ArrayList[Double]();
       coeff     = new ArrayList[Double]();
       this.type = type;

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

    public def getType() : String = this.type;
    public def getAngularMomentum() = this.angularMomentum;
    public def getExponents() : ArrayList[Double] = this.exps;
    public def getCoefficients() : ArrayList[Double] = this.coeff;
}

