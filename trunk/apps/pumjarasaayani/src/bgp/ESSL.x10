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
package org.gnu.gsl;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.Native;
import x10.compiler.NativeCPPCompilationUnit;
import x10.array.Array;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

/**
 * ESSL wrapper
 * 
 * @author V. Ganesh
 */
@NativeCPPCompilationUnit("ESSL.cc")
@NativeRep("c++", "::bgp::ESSLWrapper", "::bgp::ESSLWrapper", null)
public class ESSL {
     @Native("c++", "::bgp::ESSLWrapper::eigenSymmv((#1)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#2)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#3)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getSize())")
     public native static def eigenSymmv(A:Matrix, evec:Matrix, eval:Vector) : Int; 

     @Native("c++", "::bgp::ESSLWrapper::solve((#1)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#2)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getSize())")     
     public native static def solve(A:Matrix, b:Vector, x:Vector) : Int;
}

