/**
 * GMatrix.x10
 *
 * GMatrix in HF calculation
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

public class GMatrix extends Matrix {
    public def this(siz:Int) {
        super(siz);
    }
  
    public def compute(twoE:TwoElectronIntegrals, density:Density) : void {
        if (twoE.isDirect()) { computeDirect(twoE, density); return; }

        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector(density); // form 1D vector of density
        val tempVector  = new Vector(noOfBasisFunctions*noOfBasisFunctions);

        val gMatrix = mat;
        val ints = twoE.getTwoElectronIntegrals();
        val temp = tempVector.getVector();

        var i:Int, j:Int, k:Int, l:Int, kl:Int, 
            indexJ:Int, indexK1:Int, indexK2:Int;

        // TODO: x10 - parallel
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

    private def computeDirect(twoE:TwoElectronIntegrals, density:Density) : void {
        val noOfBasisFunctions = density.getRowCount();
        val densityOneD = new Vector(density); // form 1D vector of density
        val tempVector  = new Vector(noOfBasisFunctions*noOfBasisFunctions);

        val gMatrix = mat;
        val ints = twoE.getTwoElectronIntegrals();
        val temp = tempVector.getVector();

        var i:Int, j:Int, k:Int, l:Int, kl:Int,
            indexJ:Int, indexK1:Int, indexK2:Int;

        var twoEIntVal1:Double, twoEIntVal2:Double, twoEIntVal3:Double;

        // TODO: x10 - parallel
        for(i=0; i<noOfBasisFunctions; i++) {
            for(j=0; j<i+1; j++) {

                tempVector.makeZero();
                kl = 0;

                for(k=0; k<noOfBasisFunctions; k++) {
                    for(l=0; l<noOfBasisFunctions; l++) {
                        indexJ   = IntegralsUtils.ijkl2intindex(i, j, k, l);
                        indexK1  = IntegralsUtils.ijkl2intindex(i, k, j, l);
                        indexK2  = IntegralsUtils.ijkl2intindex(i, l, k, j);

                        twoEIntVal1 = twoE.compute2E(i,j,k,l);
                        if (indexJ == indexK1) twoEIntVal2 = twoEIntVal1;
                        else                   twoEIntVal2 = twoE.compute2E(i,k,j,l);

                        if (indexJ == indexK2)       twoEIntVal3 = twoEIntVal1;
                        else if (indexK1 == indexK2) twoEIntVal3 = twoEIntVal2;
                        else                         twoEIntVal3 = twoE.compute2E(i,l,k,j);

                        temp(kl) = 2.0*twoEIntVal1 - 0.5*twoEIntVal2 - 0.5*twoEIntVal3;

                        /**
                         x10.io.Console.OUT.println("\t" + kl);
                         x10.io.Console.OUT.println("\t\t <" + i + " " + j + " | " + k + " " + l + ">"); 
                         x10.io.Console.OUT.println("\t\t <" + i + " " + k + " | " + j + " " + l + ">"); 
                         x10.io.Console.OUT.println("\t\t <" + i + " " + l + " | " + k + " " + j + ">"); 
                        **/

                        kl++;
                    } // end l loop
                } // end k loop

                gMatrix(i, j) = gMatrix(j, i) = tempVector.dot(densityOneD);
            } // end j loop
        } // end i loop
    }
}

