/**
 * DIISFockExtrapolator.x10
 *
 * DIIS based fock extrapolation
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.matrix.Matrix;
import x10x.vector.Vector;

public class DIISFockExtrapolator {
    global val fockMatrixList:ArrayList[Fock!]!;
    global val errorMatrixList:ArrayList[Matrix!]!;

    var diisStep:Int = 0;

    public def this() {
        fockMatrixList  = new ArrayList[Fock!]();
        errorMatrixList = new ArrayList[Matrix!]();

        diisStep = 0;
    }

    public def next(currentFock:Fock!, overlap:Overlap!, density:Density!) : Fock {
        val N = currentFock.getRowCount();
        val newFock = new Fock(N) as Fock!;

        val newFockMat = newFock.getMatrix();
        val curFockMat = currentFock.getMatrix();

        for((i,j) in curFockMat.dist) newFockMat(i,j) = curFockMat(i,j);
        
        // TODO:

        val FPS = currentFock.mul(density).mul(overlap);
        val SPF = density.mul(overlap).mul(currentFock);

        val errorMatrix = FPS.sub(SPF);

        var i:Int, j:Int;

        if (diisStep > 0) {
            errorMatrixList.add(errorMatrix);
            fockMatrixList.add(currentFock);

            newFock.makeZero();

            val noOfIterations = errorMatrixList.size();
            val N1 = noOfIterations + 1;

            val A = new Matrix(N1) as Matrix!;
            val B = new Vector(N1) as Vector!;

            val aMatrix = A.getMatrix();
            val bVector = B.getVector();

            // set up A x = B to be solved
            for(i=0; i<noOfIterations; i++) {
                for(j=0; j<noOfIterations; j++) {

                    /** TODO : getRowVector() is not yet defined!
                    aMatrix(i,j) = aMatrix(j,i)
                                 = errorMatrix.getRowVector(i).dot(
                                               errorMatrix.getRowVector(j));
                    **/ 

                } // end for
            } // end for

            for(i=0; i<noOfIterations; i++) {
                aMatrix(noOfIterations,i) = aMatrix(i,noOfIterations) = -1.0;
                bVector(i) = 0.0;
            } // end for

            aMatrix(noOfIterations,noOfIterations) = 0.0;
            bVector(noOfIterations) = -1.0;

            // TODO: Call Ax=B solver, depending on the output, set the 
            // values for new Fock matrix
        } // end if

        diisStep++;

        return newFock;
    }
}

