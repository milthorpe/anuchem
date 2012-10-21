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

/**
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */
public class Density extends DenseMatrix{self.M==self.N} {
    private val noOfOccupancies:Int;

    public def this(n:Int, noOfOccupancies:Int):Density{self.M==n,self.N==n} {
        super(n, n);
        this.noOfOccupancies = noOfOccupancies;
    }

    /**
     * Creates a density matrix with each element initialised to v.
     */
    public def this(n:Int, noOfOccupancies:Int, v:Double):Density{self.M==n,self.N==n} {
        // TODO GML new constructor for DenseMatrix
        super(n, n, new Array[Double](n*n, (Int)=> v));
        this.noOfOccupancies = noOfOccupancies;
    }
/* TODO needed?
    public def this(d:Density) {
        super(d.M, d.N, d);
        this.noOfOccupancies = d.getNoOfOccupancies();
    }
*/
    public def getNoOfOccupancies() = noOfOccupancies;

    public def compute(mos:MolecularOrbitals{self.N==this.N}) : void {
        // construct it from the MOs .. C*C'
        val dVector = new DenseMatrix(noOfOccupancies, this.N);
        DenseMatrix.copyRows(mos, 0, dVector, 0, noOfOccupancies);
        super.transMult(dVector, dVector, false);
        for([x,y] in 0..(N-1)*0..(N-1)) this(x,y)=2.*this(x,y); // Closed-shell
    }

    public def applyGuess(SAD:DenseMatrix)  {
        DenseMatrix.copy(SAD, this);
    }
}

