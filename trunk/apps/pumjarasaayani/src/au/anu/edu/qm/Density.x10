/**
 * Density.x10
 *
 * The density matrix in the HF calculation
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;

public class Density extends Matrix {
    public def compute(noOfOccupancies:Int, mos:MolecularOrbitals{self.at(this)}) : void {
        // construct it from the MOs .. C*C'
        val N = mos.getRowCount();
        val dVector:Matrix{self.at(this)} = Matrix.make(noOfOccupancies, N) as Matrix{self.at(this)};

        val dMat = dVector.getMatrix();
        val mosMat = mos.getMatrix();
        for(var i:Int=0; i<noOfOccupancies; i++)
            for(var j:Int=0; j<N; j++) 
              dMat(i, j) = mosMat(i, j);

        val res = dVector.transpose().mul(dVector).getMatrix();

        val thisMat = getMatrix();
        for(val(i, j) in thisMat.region)
           thisMat(i, j) = res(i, j);
    }
}

