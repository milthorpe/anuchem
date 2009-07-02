package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

public class GMatrix extends Matrix {
    public def this(siz:Int) {
        super(siz);
    }
  
    public def compute(twoE:TwoElectronIntegrals, density:Density) : void {
        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector(density); // form 1D vector of density
        val tempVector  = new Vector(noOfBasisFunctions*noOfBasisFunctions);

        val gMatrix = mat;
        val ints = twoE.getTwoElectronIntegrals();
        val temp = tempVector.getVector();

        var i:Int, j:Int, k:Int, l:Int, kl:Int, 
            indexJ:Int, indexK1:Int, indexK2:Int;

        for(i=0; i<noOfBasisFunctions; i++) {
            for(j=0; j<i+1; j++) {

                tempVector.makeZero();
                kl = 0;

                for(k=0; k<noOfBasisFunctions; k++) {
                    for(l=0; l<noOfBasisFunctions; l++) {
                        indexJ   = IntegralsUtils.ijkl2intindex(i, j, k, l);
                        indexK1  = IntegralsUtils.ijkl2intindex(i, k, j, l);
                        indexK2  = IntegralsUtils.ijkl2intindex(i, l, k, j);
                        temp(kl) = 2.0*ints(indexJ) - 0.5*ints(indexK1)
                                   - 0.5*ints(indexK2);
                        kl++;
                    } // end l loop
                } // end k loop

                gMatrix(i, j) = gMatrix(j, i) = tempVector.dot(densityOneD);
            } // end j loop
        } // end i loop  
    }
}

