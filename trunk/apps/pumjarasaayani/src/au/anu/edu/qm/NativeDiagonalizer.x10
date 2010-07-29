/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

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
    var eigenValuesVec:Vector!;
    var eigenVectorsMat:Matrix!;

    public def diagonalize(mat:Matrix!) : void {
         val n:Int = mat.getRowCount();
         eigenVectorsMat = new Matrix(n) as Matrix!;
         eigenValuesVec  = new Vector(n) as Vector!;

         eigenVectorsMat.makeZero();
         eigenValuesVec.makeZero();

         GSL.eigenSymmv(mat, eigenVectorsMat, eigenValuesVec);
         
         eigenVectorsMat = eigenVectorsMat.transpose();
    }

    public def getEigenValues() : Vector! = eigenValuesVec;
    public def getEigenVectors() : Matrix! = eigenVectorsMat;
}

