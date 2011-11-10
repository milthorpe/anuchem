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
package au.edu.anu.qm;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.Diagonalizer;

import org.gnu.gsl.GSL;

/**
 * Native diagonalizer interface wrapper
 * 
 * @author V. Ganesh
 */
public class NativeDiagonalizer implements Diagonalizer {
    var eigenValuesVec:Vector;
    var eigenVectorsMat:Matrix;

    public def diagonalize(mat:Matrix) : void {
         val n:Int = mat.getRowCount();
         eigenVectorsMat = new Matrix(n);
         eigenValuesVec  = new Vector(n);

         GSL.eigenSymmv(mat, eigenVectorsMat, eigenValuesVec);
         
         eigenVectorsMat = eigenVectorsMat.transpose();
    }

    public def getEigenValues() : Vector = eigenValuesVec;
    public def getEigenVectors() : Matrix = eigenVectorsMat;
}

