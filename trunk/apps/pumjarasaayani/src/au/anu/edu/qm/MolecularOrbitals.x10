package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.xla.JacobiDiagonalizer;

public class MolecularOrbitals extends Matrix {
    var orbitalEnergies:Array[Double]{rank==1};

    public def this(siz:Int) {
       super(siz);
    }

    public def getOrbitalEnergies() : Array[Double]{rank==1} = orbitalEnergies;
 
    public def compute(theMat:Matrix, overlap:Overlap) : void {
        val x = overlap.symmetricOrthogonalization();
        val a = theMat.similarityTransform(x);
        val diag = new JacobiDiagonalizer();

        diag.diagonalize(a);
        orbitalEnergies = diag.getEigenValues();
        mat = diag.getEigenVectors().mul(x).getMatrix(); 
    }
}

