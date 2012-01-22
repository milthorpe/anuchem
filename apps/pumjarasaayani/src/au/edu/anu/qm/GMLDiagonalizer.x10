/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.qm;

import x10.matrix.DenseMatrix;
import x10.matrix.Vector;
import x10.matrix.lapack.DenseMatrixLAPACK;

/**
 * GML-based diagonalizer
 * 
 * @author milthorpe
 */
public class GMLDiagonalizer {
    var eigenValuesVec:Vector;
    var eigenVectorsMat:DenseMatrix;

    public def diagonalize(A:DenseMatrix):void {
        eigenValuesVec  = Vector.make(A.N) as Vector(A.N);
        val scratch = new Rail[Double](3*A.N-1);

        DenseMatrixLAPACK.compEigenVector(A, eigenValuesVec, scratch);

        // TODO required?
        eigenVectorsMat = new DenseMatrix(A.N, A.M);
        eigenVectorsMat.T(A);
    }

    public static def symmetricOrthogonalization(A:DenseMatrix):DenseMatrix(A.N,A.N) {
        val diag = new GMLDiagonalizer();
        diag.diagonalize(A);
        
        val eigenValues   = diag.getEigenValues();
        val eigenVectors  = diag.getEigenVectors();
        val sHalf         = new DenseMatrix(A.M, A.N);

        // TODO a smarter way to do this...
        for (i in 0..(eigenValues.N-1)) {
            sHalf(i,i) = 1.0 / Math.sqrt(eigenValues(i));
        }

        val x = new DenseMatrix(A.M, A.N);
        x.transMult(eigenVectors, sHalf % eigenVectors);
        return x;
    }

    public def getEigenValues():Vector = eigenValuesVec;
    public def getEigenVectors():DenseMatrix = eigenVectorsMat;
}

