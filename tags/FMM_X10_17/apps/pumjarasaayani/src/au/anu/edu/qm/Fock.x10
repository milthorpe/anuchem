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
    public def compute(hCore:HCore, gMatrix:GMatrix) : void {
        val res = hCore.add(gMatrix).getMatrix();
        val thisMat = getMatrix();

        for(val(i, j) in res.region)
           thisMat(i, j) = res(i, j);
    }
}

