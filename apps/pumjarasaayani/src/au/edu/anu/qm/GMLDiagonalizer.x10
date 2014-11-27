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
    var eigenvaluesVec:Vector;
    var eigenvectorsMat:DenseMatrix{self.M==self.N};

    public def diagonalize(A:DenseMatrix{self.M==self.N}):void {
        eigenvectorsMat = new DenseMatrix(A.N,A.M);
        eigenvaluesVec = Vector.make(eigenvectorsMat.N);

        DenseMatrixLAPACK.compEigenvectors(A.clone(), eigenvaluesVec, eigenvectorsMat);
    }

    public static def symmetricOrthogonalization(A:DenseMatrix{self.M==self.N}):DenseMatrix(A.N,A.N) {
        val diag = new GMLDiagonalizer();
        diag.diagonalize(A);

        val eigenvalues   = diag.getEigenvalues();
        val eigenvectors  = diag.getEigenvectors() as DenseMatrix(A.N,A.M);
        val sHalf         = new DenseMatrix(A.M, A.N);

        // TODO a smarter way to do this...
        for (i in 0..(eigenvalues.M-1)) {
            sHalf(i,i) = 1.0 / Math.sqrt(eigenvalues(i));
        }

        val x = new DenseMatrix(A.M, A.N);
        x.multTrans(eigenvectors % sHalf, eigenvectors);
        return x;
    }

    public def getEigenvalues():Vector = eigenvaluesVec;
    public def getEigenvectors():DenseMatrix = eigenvectorsMat;
}

