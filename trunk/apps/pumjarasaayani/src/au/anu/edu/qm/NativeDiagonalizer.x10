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
    global var eigenValuesVec:Vector{self.at(this)};
    global var eigenVectorsMat:Matrix{self.at(this)};

    public def diagonalize(mat:Matrix{self.at(this)}) : void {
         val n:Int = mat.getRowCount();
         eigenVectorsMat = Matrix.make(n) as Matrix{self.at(this)};
         eigenValuesVec  = Vector.make(n) as Vector{self.at(this)};

         eigenVectorsMat.makeZero();
         eigenValuesVec.makeZero();

         GSL.eigenSymmv(mat, eigenVectorsMat, eigenValuesVec);
         
         eigenVectorsMat = eigenVectorsMat.transpose();
    }

    public def getEigenValues() : Vector{self.at(this)} = eigenValuesVec;
    public def getEigenVectors() : Matrix{self.at(this)} = eigenVectorsMat;
}

