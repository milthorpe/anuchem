/**
 * MolecularOrbitals.x10
 *
 * Represents MOs in a HF-SCF
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.xla.JacobiDiagonalizer;

public class MolecularOrbitals extends Matrix {
    var orbitalEnergies:Array[Double]{rank==1};

    public def getOrbitalEnergies() : Array[Double]{rank==1} = orbitalEnergies;
 
    public def compute(theMat:Matrix{self.at(this)}, overlap:Overlap{self.at(this)}) : void {
        val x:Matrix{self.at(this)} = overlap.getSHalf();
        val a:Matrix{self.at(this)} = theMat.similarityTransform(x);
        val diag = new JacobiDiagonalizer();
        // val diag = new NativeDiagonalizer();

        diag.diagonalize(a);
        orbitalEnergies = diag.getEigenValues().getVector();
        val res = diag.getEigenVectors().mul(x).getMatrix(); 
        val thisMat = getMatrix();

        for(val(i, j) in res.region)
           thisMat(i, j) = res(i, j);
    }
}

