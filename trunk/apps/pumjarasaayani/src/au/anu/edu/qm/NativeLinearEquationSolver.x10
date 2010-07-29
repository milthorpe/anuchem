/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

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

     public def this() {}
 
     public def findSolution(matrixA:Matrix!, vectorB:Vector!) : Vector! throws Exception {
          val N = matrixA.getRowCount();
          val vectorX = new Vector(N) as Vector!;

          vectorX.makeZero();
          
          if (GSL.solve(matrixA, vectorB, vectorX) != 0) throw new Exception("No solution!");

          return vectorX;
     }
}

