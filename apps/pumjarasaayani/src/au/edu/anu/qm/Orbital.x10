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

/**
 * Represents an Orbital (used for basis function formation)
 *
 * @author: V.Ganesh, milthorpe
 */
public struct Orbital {
    public val exponents:Rail[Double];
    public val coefficients:Rail[Double];
    public val shape:String;
    public val angularMomentum:Int;

    public def this(shape:String, exps:Rail[Double], coeffs:Rail[Double]) { 
       this.exponents = exps;
       this.coefficients = coeffs;
       this.shape = shape;

       if (shape.equals("S")) angularMomentum = 0n;
       else if (shape.equals("P")) angularMomentum = 1n;
       else if (shape.equals("D")) angularMomentum = 2n;
       else if (shape.equals("F")) angularMomentum = 3n;
       else if (shape.equals("G")) angularMomentum = 4n;
       else if (shape.equals("H")) angularMomentum = 5n;
//     else if (shape.equals("I")) angularMomentum = 6n;
       else {Console.OUT.printf("Orbital.x10 support SPDFGH: unknown orbital '%s'.\n",shape); angularMomentum = 0n;}
    }
}

