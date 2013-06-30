/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.util.ArrayList;
import x10.matrix.DenseMatrix;
import x10.matrix.Vector;
import x10.matrix.lapack.DenseMatrixLAPACK;

/**
 * DIIS based fock extrapolation
 * 
 * @see Pulay, P. (1980). "Convergence acceleration of iterative sequences.
 *     The case of SCF iteration".  Chem. Phys. Lett., 73 (2) pp. 393-398.
 * @see Pulay, P. (1982). "Improved SCF convergence acceleration".
 *     J. Comp. Chem., 3 (4) pp. 556-560.
 *
 * @author: V.Ganesh
 */
public class DIISFockExtrapolator(N:Long) {
    static MIN_NON_DIIS_STEP:Int = 0;
    static MAX_NON_DIIS_STEP:Int = 2;

    private val diisStartThreshold:Double;
    private val diisSubspaceSize:Int;
    private val diisMaxError:Double;

    val fockMatrixList:ArrayList[Fock];
    val errorVectorList:ArrayList[Vector];
    val FPS:DenseMatrix{self.M==self.N,self.N==this.N};
    val SPF:DenseMatrix{self.M==self.N,self.N==this.N};
    val scratch:DenseMatrix{self.M==self.N,self.N==this.N};

    var diisStep:Int = 0;
    var nondiisStep:Int = 0;

    var diisStarted:Boolean = false;
    var converged:Boolean = false;

    var oldFock:Fock;

    public def this(N:Long) {
        property(N);
        fockMatrixList  = new ArrayList[Fock]();
        errorVectorList = new ArrayList[Vector]();

        FPS = new DenseMatrix(N,N);
        SPF = new DenseMatrix(N,N);
        scratch = new DenseMatrix(N,N);

        val jd = JobDefaults.getInstance();
        diisStartThreshold = jd.diisStartThreshold;
        diisMaxError = jd.diisConvergenceThreshold;
        diisSubspaceSize = jd.diisSubspaceSize;

        nondiisStep = 0;
        diisStep = 0;
        converged = false;
    }

    public def isConverged() : Boolean {
        return converged;
    }

    /** Generate a new Fock by extrapolating from previous recorded Fock and their difference vectors */
    public def next(currentFock:Fock{self.M==this.N,self.N==this.N}, 
                    overlap:Overlap{self.M==currentFock.M,self.N==currentFock.N},
                    density:Density{self.M==currentFock.M,self.N==currentFock.N}):Fock{self.M==self.N,self.N==currentFock.N} {
        val newFock = new Fock(N);

        scratch.mult(currentFock, density);
        FPS.mult(scratch, overlap);
        scratch.mult(overlap, density);
        SPF.mult(scratch, currentFock);

        val errorVector = new Vector((FPS - SPF).d);
        val mxerr = errorVector.maxNorm();
        val errorVectorSize = errorVectorList.size();

        Console.OUT.printf("   Max DIIS error = %.3e subspaceSize = %d\n",mxerr,errorVectorSize);
        if (mxerr < diisMaxError) {
            converged=true;
            return currentFock;
        }

        if (errorVectorSize >= diisSubspaceSize) {
            errorVectorList.removeAt(0);
            fockMatrixList.removeAt(0);
        }

        // TODO if (reset=true) remove all vector
        
        errorVectorList.add(errorVector);
        fockMatrixList.add(currentFock);

        val noOfIterations = errorVectorList.size();

        if (!diisStarted && 
          ( (mxerr < diisStartThreshold && diisStep>=MIN_NON_DIIS_STEP) || nondiisStep>=MAX_NON_DIIS_STEP)) {
            Console.OUT.println("Starting DIIS...");
            diisStarted = true;
        }

        if (!diisStarted) {
            nondiisStep++;
            if (oldFock == null) {
                oldFock = currentFock;

                return currentFock;
            } else {
                for(var i:Int=0; i<N; i++) {
                   for(var j:Int=0; j<N; j++) {
                      newFock(i, j) = oldFock(i,j) * 0.5 + currentFock(i,j) * 0.5;
                   }
                }

                oldFock = currentFock;

                return newFock;
            } // end if
        } // end if 

        val M = noOfIterations+1;
        val A = new DenseMatrix(M,M);
        val B = Vector.make(M);

        // set up A x = B to be solved
        for (var i:Int=0; i < noOfIterations; i++) {
            for (var j:Int=0; j < noOfIterations; j++) {
                A(i,j) = errorVectorList.get(i).blasTransProduct(errorVectorList.get(j));
            } // end for
        } // end for

        for (var i:Int=0; i < noOfIterations; i++) {
            A(noOfIterations,i) = A(i,noOfIterations) = -1.0;
            B(i) = 0.0;
        } // end for

        A(noOfIterations,noOfIterations) = 0.0;
        B(noOfIterations) = -1.0;

        val permutation = new Rail[Int](M);
        val result = DenseMatrixLAPACK.solveLinearEquation(A, B, permutation);

        for (var i:Int=0; i < noOfIterations; i++) {
          val prevFock = fockMatrixList.get(i);
          for (var j:Int=0; j < N; j++) {
             for (var k:Int=0; k < N; k++) {
                 newFock(j,k) += B(i) * prevFock(j,k);
             } // end for
          } // end for
        } // end for

        diisStep++;

        return newFock;
    }
}

