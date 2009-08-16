/**
 * Fock.x10
 *
 * The Fock matrix in the HF calculation
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;

public class Fock extends Matrix {
    public def this(siz:Int) {
       super(siz);
    }

    public def compute(hCore:HCore, gMatrix:GMatrix) : void {
        val res = hCore.add(gMatrix).getMatrix();
        val N   = getRowCount();
        
        var i:Int, j:Int;
        for(i=0; i<N; i++)
           for(j=0; j<N; j++)
              mat(i, j) = res(i, j);
    }
}

