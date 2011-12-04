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

import x10x.matrix.Matrix;

/**
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Density extends Matrix {
    private val noOfOccupancies:Int;

    public def this(n:Int, noOfOccupancies:Int) {
        super(n);
        this.noOfOccupancies = noOfOccupancies;
    }

    public def this(d:Density) {
        super(d);
        this.noOfOccupancies = d.getNoOfOccupancies();
    }

    public def getNoOfOccupancies() = noOfOccupancies;

    public def compute(mos:MolecularOrbitals) : void {
        // construct it from the MOs .. C*C'
        val N = mos.getRowCount();
        val dVector = new Matrix(noOfOccupancies, N);

        val dMat = dVector.getMatrix();
        val mosMat = mos.getMatrix();
        for(var i:Int=0; i<noOfOccupancies; i++)
            for(var j:Int=0; j<N; j++) 
              dMat(i, j) = mosMat(i, j);

        this.mulInPlace(dVector.transpose(), dVector);
    }

    public def applyGuess(SAD:Matrix)  {
        val N = SAD.getRowCount();
        val dMat = SAD.getMatrix();
        val thisMat = getMatrix();
        Array.copy[Double](dMat, thisMat);
    }
}

