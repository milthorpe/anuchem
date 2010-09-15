/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.anu.edu.qm;

import x10.util.ArrayList;
import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.GaussianElimination;

/**
 * DIISFockExtrapolator.x10
 *
 * DIIS based fock extrapolation
 * Note: Mostly lifted from MeTA Studio code.
 * For reference see:
 * "CONVERGENCE ACCELERATION OF ITERATIVE SEQUENCES. THE CASE OF SCF ITERATION", Peter Pulay, Chem. Phys. Lett., 73, 393, 1980.
 * (opps! thats my birth year ;-))
 *
 * @author: V.Ganesh
 */
public class DIISFockExtrapolator {
    val fockMatrixList:ArrayList[Fock];
    val errorVectorList:ArrayList[Vector];

    val errorThreshold = 0.1;

    var diisStep:Int = 0;
    var diisStarted:Boolean = false;

    var oldFock:Fock;

    public def this() {
        fockMatrixList  = new ArrayList[Fock]();
        errorVectorList = new ArrayList[Vector]();

        diisStep = 0;
    }

    /** Generate a new Fock by extrapolating from previous recorded Fock and their difference vectors */
    public def next(currentFock:Fock, overlap:Overlap, density:Density) : Fock {
        val N = currentFock.getRowCount();
        var newFock:Fock = new Fock(N);

        val newFockMat = newFock.getMatrix();
        val curFockMat = currentFock.getMatrix();

        var i:Int, j:Int, k:Int;

        for(i=0; i<N; i++)
           for(j=0; j<N; j++)
              newFockMat(i,j) = curFockMat(i,j);
        
        val FPS = currentFock.mul(density).mul(overlap);
        val SPF = overlap.mul(density).mul(currentFock);

        val errorVector = new Vector(FPS.sub(SPF));
        val mxerr = errorVector.maxNorm();

        if (mxerr < errorThreshold && !diisStarted) {
            Console.OUT.println("Starting DIIS...");
            diisStarted = true;
        } // end if

        if (!diisStarted) {
            if (oldFock == null) {
                oldFock = currentFock;

                return currentFock;
            } else {
                val nf = oldFock.mul(0.5).add(currentFock.mul(0.5)).getMatrix();
                for(i=0; i<N; i++) {
                   for(j=0; j<N; j++) {
                      newFockMat(i, j) = nf(i, j);
                   }
                }

                oldFock = currentFock;

                return newFock;
            } // end if
        } // end if
        
        errorVectorList.add(errorVector);
        fockMatrixList.add(currentFock);

        newFock.makeZero();

        val noOfIterations = errorVectorList.size();

        val A = new Matrix(noOfIterations+1);
        val B = new Vector(noOfIterations+1);

        val aMatrix = A.getMatrix();
        val bVector = B.getVector();

        // set up A x = B to be solved
        for (i = 0; i < noOfIterations; i++) {
            for (j = 0; j < noOfIterations; j++) {
                aMatrix(i,j) = errorVectorList.get(i).dot(errorVectorList.get(j));
            } // end for
        } // end for

        for (i = 0; i < noOfIterations; i++) {
            aMatrix(noOfIterations,i) = aMatrix(i,noOfIterations) = -1.0;
            bVector(i) = 0.0;
        } // end for

        aMatrix(noOfIterations,noOfIterations) = 0.0;
        bVector(noOfIterations) = -1.0;

        // val gele = new GaussianElimination();
        val gele = new NativeLinearEquationSolver();

        try {
          val solVec = gele.findSolution(A, B).getVector();

          for (i = 0; i < noOfIterations; i++) {
              val prevFockMat = fockMatrixList.get(i).getMatrix();
              for (j = 0; j < N; j++) {
                 for (k = 0; k < N; k++) {
                     newFockMat(j,k) += solVec(i) * prevFockMat(j,k);
                 } // end for
              } // end for
          } // end for
        } catch(e:Exception) {
          diisStep++;
          
          return currentFock;
        } // end catch throw

        diisStep++;

        return newFock;
    }
}

