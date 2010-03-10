/**
 * Native diagonalizer interface wrapper 
 *
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.Diagonalizer;

import org.gnu.gsl.GSL;

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

