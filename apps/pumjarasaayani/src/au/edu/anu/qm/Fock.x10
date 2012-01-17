/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10x.matrix.Matrix;

/**
 * The Fock matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Fock extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    public def compute(hCore:HCore, gMatrix:Matrix) : void {
        this.addInPlace(hCore, gMatrix);
    }
}

