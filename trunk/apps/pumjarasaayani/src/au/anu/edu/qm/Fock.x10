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

import x10x.matrix.Matrix;

/**
 * Fock.x10
 *
 * The Fock matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Fock extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    public def compute(hCore:HCore, gMatrix:GMatrix) : void {
        val res = hCore.add(gMatrix).getMatrix();
        val thisMat = getMatrix();

        for([i, j] in res)
           thisMat(i, j) = res(i, j);
    }
}

