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

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.LinearEquationSolver;

import org.gnu.gsl.GSL;

/**
 * Native Ax=B solver interface wrapper
 *
 * @author V. Ganesh
 */
public class NativeLinearEquationSolver extends LinearEquationSolver {
 
    public def findSolution(matrixA:Matrix, vectorB:Vector) : Vector {
        val N = matrixA.getRowCount();
        val vectorX = new Vector(N);

        val solvedResult = GSL.solve(matrixA, vectorB, vectorX);
        if (solvedResult != 0) {
            Console.ERR.println("NativeLinearEquationSolver.findSolution(): no solution! return code = " + solvedResult);
        }

        return vectorX;
    }
}

