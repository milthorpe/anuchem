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
import x10x.matrix.Matrix;

/**
 * Density.x10
 *
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Density extends Matrix {
    public def this(n:Int) {
        super(n);
    }

    public def compute(noOfOccupancies:Int, mos:MolecularOrbitals!) : void {
        // unmute();

        // construct it from the MOs .. C*C'
        val N = mos.getRowCount();
        val dVector:Matrix! = new Matrix(noOfOccupancies, N) as Matrix!;

        val dMat = dVector.getMatrix();
        val mosMat = mos.getMatrix();
        for(var i:Int=0; i<noOfOccupancies; i++)
            for(var j:Int=0; j<N; j++) 
              dMat(i, j) = mosMat(i, j);

        val res = dVector.transpose().mul(dVector).getMatrix();

        val thisMat = getMatrix();
        for(val(i, j) in thisMat.region)
           thisMat(i, j) = res(i, j);

        // mute();
    }
}

