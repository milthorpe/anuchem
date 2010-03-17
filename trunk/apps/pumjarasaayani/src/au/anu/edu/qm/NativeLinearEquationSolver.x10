/**
 * Native Ax=B solver interface wrapper
 *
 */

package au.anu.edu.qm;

import x10x.matrix.Matrix;
import x10x.vector.Vector;
import x10x.xla.LinearEquationSolver;

import org.gnu.gsl.GSL;

public class NativeLinearEquationSolver extends LinearEquationSolver {

     public def this() {}
 
     public def findSolution(matrixA:Matrix!, vectorB:Vector!) : Vector! throws Exception {
          val N = matrixA.getRowCount();
          val vectorX = new Vector(N) as Vector!;

          vectorX.makeZero();
          
          if (GSL.solve(matrixA, vectorB, vectorX) != 0) throw new Exception("No solution!");

          return vectorX;
     }
}

