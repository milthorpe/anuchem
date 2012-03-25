/**
 * An abstract class representing a linear equation solver of the form: Ax = B
 *
 * @author  V.Ganesh
 */

package x10x.xla;

import x10x.vector.Vector;
import x10x.matrix.Matrix;

public abstract class LinearEquationSolver {

    public def this() {
    }

    public abstract def findSolution(matrixA:Matrix, vectorB:Vector) : Vector;
}

