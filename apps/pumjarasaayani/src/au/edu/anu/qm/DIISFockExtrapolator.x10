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
package au.edu.anu.qm;

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
    val MinnondiisStep:Int = 2;
    val MaxnondiisStep:Int = 5;    

    var diisStep:Int = 0;
    var nondiisStep:Int = 0;

    var diisStarted:Boolean = false;
    var converged:Boolean = false;

    var oldFock:Fock;

    public def this() {
        fockMatrixList  = new ArrayList[Fock]();
        errorVectorList = new ArrayList[Vector]();
        nondiisStep = 0;
        diisStep = 0;
        converged = false;
    }

    public def isConverged() : Boolean {
        return converged;
    }

    /** Generate a new Fock by extrapolating from previous recorded Fock and their difference vectors */
    public def next(currentFock:Fock, overlap:Overlap, density:Density) : Fock {
        val N = currentFock.getRowCount();
        var newFock:Fock = new Fock(N);

        val newFockMat = newFock.getMatrix();
        val curFockMat = currentFock.getMatrix();
        
        val FPS = currentFock.mul(density).mul(overlap);
        val SPF = overlap.mul(density).mul(currentFock);

        val errorVector = new Vector(FPS.sub(SPF));
        val mxerr = errorVector.maxNorm();

        //Console.OUT.printf("Max DIIS error %20.11f\n",mxerr);
        if (mxerr <1e-6) converged=true;

        if ( ( (mxerr < errorThreshold && diisStep>=MaxnondiisStep) || nondiisStep>=MaxnondiisStep) && !diisStarted) {
            Console.OUT.println("Starting DIIS...");
            diisStarted = true;
        } // end if

        if (!diisStarted) {
            nondiisStep++;
            if (oldFock == null) {
                oldFock = currentFock;

                return currentFock;
            } else {
                val oldFockMat = oldFock.getMatrix();
                for(var i:Int=0; i<N; i++) {
                   for(var j:Int=0; j<N; j++) {
                      newFockMat(i, j) = oldFockMat(i,j) * 0.5 + curFockMat(i,j) * 0.5;
                   }
                }

                oldFock = currentFock;

                return newFock;
            } // end if
        } // end if
        
        errorVectorList.add(errorVector);
        fockMatrixList.add(currentFock);

        val noOfIterations = errorVectorList.size();

        val A = new Matrix(noOfIterations+1);
        val B = new Vector(noOfIterations+1);

        val aMatrix = A.getMatrix();
        val bVector = B.getVector();

        // set up A x = B to be solved
        for (var i:Int=0; i < noOfIterations; i++) {
            for (var j:Int=0; j < noOfIterations; j++) {
                aMatrix(i,j) = errorVectorList.get(i).dot(errorVectorList.get(j));
            } // end for
        } // end for

        for (var i:Int=0; i < noOfIterations; i++) {
            aMatrix(noOfIterations,i) = aMatrix(i,noOfIterations) = -1.0;
            bVector(i) = 0.0;
        } // end for

        aMatrix(noOfIterations,noOfIterations) = 0.0;
        bVector(noOfIterations) = -1.0;

        // val gele = new GaussianElimination();
        val gele = new NativeLinearEquationSolver();

        val solVec = gele.findSolution(A, B).getVector();

        for (var i:Int=0; i < noOfIterations; i++) {
          val prevFockMat = fockMatrixList.get(i).getMatrix();
          for (var j:Int=0; j < N; j++) {
             for (var k:Int=0; k < N; k++) {
                 newFockMat(j,k) += solVec(i) * prevFockMat(j,k);
             } // end for
          } // end for
        } // end for

        diisStep++;

        return newFock;
    }
}

