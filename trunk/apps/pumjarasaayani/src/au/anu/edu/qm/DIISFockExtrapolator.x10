/**
 * DIISFockExtrapolator.x10
 *
 * DIIS based fock extrapolation
 * Note: Mostly lifted from MeTA Studio code, except no Matrix singularity checks made!
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.GaussianElimination;

public class DIISFockExtrapolator {
    global val fockMatrixList:ArrayList[Fock!]!;
    global val errorMatrixList:ArrayList[Vector!]!;

    global val errorThreshold = 0.1;

    var diisStep:Int = 0;
    var diisStarted:Boolean = false;

    var oldFock:Fock!;

    public def this() {
        fockMatrixList  = new ArrayList[Fock!]();
        errorMatrixList = new ArrayList[Vector!]();

        diisStep = 0;
    }

    public def next(currentFock:Fock!, overlap:Overlap!, density:Density!) : Fock! {
        val N = currentFock.getRowCount();
        val newFock = new Fock(N) as Fock!;

        var newFockMat:Array[Double]{rank==2, self.at(this)} = newFock.getMatrix() as Array[Double]{rank==2, self.at(this)};
        val curFockMat = currentFock.getMatrix();

        for((i,j) in curFockMat.dist) newFockMat(i,j) = curFockMat(i,j);
        
        val FPS = currentFock.mul(density).mul(overlap);
        val SPF = density.mul(overlap).mul(currentFock);

        val errorMatrix = new Vector(FPS.sub(SPF)) as Vector!;
        val mxerr = errorMatrix.maxNorm();

        if (mxerr < errorThreshold && !diisStarted) {
            Console.OUT.println("Starting DIIS...");
            diisStarted = true;
        } // end if

        if (!diisStarted) {
            if (oldFock == null) {
                oldFock = currentFock;

                return currentFock;
            } else {
                newFockMat = oldFock.mul(0.5).add(currentFock.mul(0.5)).getMatrix() as Array[Double]{rank==2, self.at(this)};
                oldFock = currentFock;

                return newFock;
            } // end if
        } // end if
        
        errorMatrixList.add(errorMatrix);
        fockMatrixList.add(currentFock);

        newFock.makeZero();

        val noOfIterations = errorMatrixList.size();
        val N1 = noOfIterations + 1;

        val A = new Matrix(N1);
        val B = new Vector(N1);

        val aMatrix = A.getMatrix();
        val bVector = B.getVector();

        var i:Int, j:Int, k:Int;

        // set up A x = B to be solved
        for (i = 0; i < noOfIterations; i++) {
            for (j = 0; j < noOfIterations; j++) {
                aMatrix(i,j) = errorMatrixList.get(i).dot(errorMatrixList.get(j));
            } // end for
        } // end for

        for (i = 0; i < noOfIterations; i++) {
            aMatrix(noOfIterations,i) = aMatrix(i,noOfIterations) = -1.0;
            bVector(i) = 0.0;
        } // end for

        aMatrix(noOfIterations,noOfIterations) = 0.0;
        bVector(noOfIterations) = -1.0;

        val gele = new GaussianElimination();

        val solVec = gele.findSolution(A, B).getVector();

        for (i = 0; i < noOfIterations; i++) {
            val prevFockMat = fockMatrixList.get(i).getMatrix();
            for (j = 0; j < N; j++) {
                 for (k = 0; k < N; k++) {
                     newFockMat(j,k) += solVec(i) * prevFockMat(j,k);
                 } // end for
            } // end for
        } // end for

        diisStep++;

        return newFock;
    }
}

