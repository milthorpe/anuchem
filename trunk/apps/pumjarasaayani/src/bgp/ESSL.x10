package org.gnu.gsl;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.Native;
import x10.compiler.NativeCPPCompilationUnit;
import x10.array.Array;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

@NativeCPPCompilationUnit("ESSL.cc")
@NativeRep("c++", "::bgp::ESSLWrapper", "::bgp::ESSLWrapper", null)
public class ESSL {
     @Native("c++", "::bgp::ESSLWrapper::eigenSymmv((#1)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#2)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#3)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getSize())")
     public native static def eigenSymmv(A:Matrix, evec:Matrix, eval:Vector) : Int; 

     @Native("c++", "::bgp::ESSLWrapper::solve((#1)._val->getMatrix()._val->raw()->FMGL(chunk)->data, (#2)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getVector()._val->raw()->FMGL(chunk)->data, (#3)._val->getSize())")     
     public native static def solve(A:Matrix, b:Vector, x:Vector) : Int;
}

