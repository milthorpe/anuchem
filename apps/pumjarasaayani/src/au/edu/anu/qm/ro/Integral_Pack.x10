/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2012.
 */
package au.edu.anu.qm.ro;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.NativeCPPCompilationUnit;
import x10.compiler.NativeCPPInclude;
import x10.array.Array;

/**
 * Wrapper for C code for auxiliary integral calculation
 *
 * @author milthorpe
 */
@NativeCPPCompilationUnit("bessel4.cc")
@NativeCPPCompilationUnit("Integral_Pack.cc")
@NativeRep("c++", "::au::edu::anu::qm::ro::Integral_Pack *", "::au::edu::anu::qm::ro::Integral_Pack", null)
public class Integral_Pack {

    public native def this(N:Int,L:Int);

    /**
     * Generate a class of auxiliary integrals for the given parameters
     */
    @Native("c++", "(#this)->Genclass((#1), (#2), (#3)._val->raw().raw(), (#4)._val->raw().raw(), (#5)._val->raw().raw(), (#6)._val->raw().raw(), (#7)._val->raw().raw(), (#8)._val->raw().raw(), (#9), (#10), (#11)._val->raw().raw())")
    public native def genClass(a:Int, b:Int, A:Rail[Double], B:Rail[Double], zetaA:Rail[Double], zetaB:Rail[Double], conA:Rail[Double], conB:Rail[Double], dconA:Int, dconB:Int, temp:Array[Double]/*(2)*/{rect,zeroBased}):Int;

}

