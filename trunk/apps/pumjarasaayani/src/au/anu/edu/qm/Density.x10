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
    public def this(siz:Int) {
       super(siz);
    }

    public def compute(noOfOccupancies:Int, mos:MolecularOrbitals) : void {
        // construct it from the MOs .. C*C'
        val N = mos.getRowCount();
        val dVector = new Matrix(noOfOccupancies, N);
        var i:Int, j:Int;

        // TODO : x10 - parallel
        for(i=0; i<noOfOccupancies; i++)
            for(j=0; j<N; j++)
                dVector.mat(i, j) = mos.mat(i, j);

        val res = dVector.transpose().mul(dVector).getMatrix();
        for(i=0; i<N; i++)
           for(j=0; j<N; j++)
              mat(i, j) = res(i, j);
    }
}

