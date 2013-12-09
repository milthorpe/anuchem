/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package org.gnu.gsl;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.Native;
import x10.compiler.NativeCPPCompilationUnit;
import x10.array.Array;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

/**
 * Wrapper for GSL
 *
 * @author V. Ganesh
 */
@NativeCPPCompilationUnit("GSL.cc")
@NativeRep("c++", "org::gnu::gsl::GSLWrapper", "org::gnu::gsl::GSLWrapper", null)
public class GSL {
    /**
     * Computes the eigenvalues and eigenvectors of a real symmetric matrix.
     * @param A a real symmetric matrix
     * @param evec a Matrix of the same dimensions as A; on return each row of 
     *   evec will contain an eigenvector of A corresponding to an eigenvalue
     * @param eval a vector containing the eigenvalues of A in ascending order
     */
    @Native("c++", "::org::gnu::gsl::GSLWrapper::eigenSymmv((#1)._val->getMatrix()._val->raw()->raw(), (#2)._val->getMatrix()._val->raw()->raw(), (#3)._val->getVector()._val->raw()->raw(), (#3)._val->getSize())")
    public native static def eigenSymmv(A:Matrix, evec:Matrix, eval:Vector) : Int; 

    /**
     * Solve the system of linear equations Ax = b using LU decomposition.
     * @param A the matrix A
     * @param b the result vector B
     * @param x the unknown vector x which is to be solved
     */
    @Native("c++", "::org::gnu::gsl::GSLWrapper::solve((#1)._val->getMatrix()._val->raw()->raw(), (#2)._val->getVector()._val->raw()->raw(), (#3)._val->getVector()._val->raw()->raw(), (#3)._val->getSize())")     
    public native static def solve(A:Matrix, b:Vector, x:Vector) : Int;
}

