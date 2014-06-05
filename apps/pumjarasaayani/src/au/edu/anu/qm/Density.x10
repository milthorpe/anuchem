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

import x10.matrix.DenseMatrix;
import x10.matrix.blas.DenseMatrixBLAS;

/**
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Density extends DenseMatrix{self.M==self.N} {
    private val noOfOccupancies:Long;

    public def this(n:Long, noOfOccupancies:Long):Density{self.M==n,self.N==n} {
        super(n, n);
        this.noOfOccupancies = noOfOccupancies;
    }

    /**
     * Creates a density matrix with each element initialised to v.
     */
    public def this(n:Long, noOfOccupancies:Long, v:Double):Density{self.M==n,self.N==n} {
        // TODO GML new constructor for DenseMatrix
        super(n, n, new Rail[Double](n*n, v));
        this.noOfOccupancies = noOfOccupancies;
    }

    public def getNoOfOccupancies() = noOfOccupancies;

    public def compute(mos:MolecularOrbitals{self.N==this.N}) : void {
        // construct it from the MOs .. C*C'
        val offsets = new Rail[Long](6);
        DenseMatrixBLAS.compTransMult(2.0, mos, mos, 0.0, this, [mos.N, mos.N, noOfOccupancies], offsets);
    }

    public def applyGuess(SAD:DenseMatrix)  {
        DenseMatrix.copy(SAD, this);
    }
}

