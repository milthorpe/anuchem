/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2013.
 */
package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

/**
 * An abstract class representing a linear equation solver of the form: Ax = B
 *
 * @author  V.Ganesh
 */
public abstract class LinearEquationSolver {

    public abstract def findSolution(matrixA:Matrix, vectorB:Vector) : Vector;
}

