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

import x10.util.*;

/**
 * Orbitals.x10
 *
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

       if (shape.equals("S")) angularMomentum = 0;
       else if (shape.equals("P")) angularMomentum = 1;
       else if (shape.equals("D")) angularMomentum = 2;
       else if (shape.equals("F")) angularMomentum = 3;
       else angularMomentum = 0;
    }
}

