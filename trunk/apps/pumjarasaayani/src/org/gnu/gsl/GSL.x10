package org.gnu.gsl;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.array.BaseArray;

import x10x.matrix.Matrix;
import x10x.vector.Vector;

@NativeRep("c++", "::org::gnu::gsl::GSLWrapper", "::org::gnu::gsl::GSLWrapper", null)
public class GSL {
     @Native("c++", "::org::gnu::gsl::GSLWrapper::eigenSymmv((reinterpret_cast<x10::array::BaseArray<double>*>(#1._val->getMatrix()._val))->raw()._val->_data, (reinterpret_cast<x10::array::BaseArray<double>*>(#2._val->getMatrix()._val))->raw()._val->_data, (reinterpret_cast<x10::array::BaseArray<double>*>(#3._val->getVector()._val))->raw()._val->_data, #3._val->getSize())")
     public native static def eigenSymmv(A:Matrix, evec:Matrix, eval:Vector) : Int; 

     @Native("c++", "::org::gnu::gsl::GSLWrapper::solve((reinterpret_cast<x10::array::BaseArray<double>*>(#1._val->getMatrix()._val))->raw()._val->_data, (reinterpret_cast<x10::array::BaseArray<double>*>(#2._val->getVector()._val))->raw()._val->_data, (reinterpret_cast<x10::array::BaseArray<double>*>(#3._val->getVector()._val))->raw()._val->_data, #3._val->getSize())")     
     public native static def solve(A:Matrix, b:Vector, x:Vector) : Int;
}

