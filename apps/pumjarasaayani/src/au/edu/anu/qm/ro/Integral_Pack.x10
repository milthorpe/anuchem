/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2012-2013.
 */
package au.edu.anu.qm.ro;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.NativeCPPCompilationUnit;
import x10.compiler.NativeCPPInclude;

import x10x.vector.Point3d;

/**
 * Wrapper for C code for auxiliary integral calculation
 *
 * @author milthorpe
 */
@NativeRep("c++", "au::edu::anu::qm::ro::Integral_Pack *", "au::edu::anu::qm::ro::Integral_Pack", null)
public class Integral_Pack {

    //@Native("c++","(#this)->Integral_Pack->_make((#1),(#2),(#3),(#4),(#5),(#6)->raw)") // doesn't work
    public native def this(N:Int,L:Int,Type:Double,roThresh:Double,rad:Double,roZ:Double);

    /**
     * Generate Ylm for a given shell pair
     */
    @Native("c++", "(#this)->GenclassY((double*)&((#1).x10__i), (double*)&((#2).x10__i), (#3)->raw, (#4)->raw, (#5), (#6), (#7), (#8)->raw)")
    public native def genClassY(A:Point3d, B:Point3d, zetaA:Rail[Double], zetaB:Rail[Double], dconA:Int, dconB:Int, Ln:Int, Ylm:Rail[Double]):Int;

    /**
     * Generate a class of auxiliary integrals for the given parameters
     */
    @Native("c++", "(#this)->Genclass((#1), (#2), (double*)&((#3).x10__i), (double*)&((#4).x10__i), (#5)->raw, (#6)->raw, (#7)->raw, (#8)->raw, (#9), (#10), (#11)->raw, (#12), (#13), (#14)->raw, (#15))")
    public native def genClass(a:Int, b:Int, A:Point3d, B:Point3d, zetaA:Rail[Double], zetaB:Rail[Double], conA:Rail[Double], conB:Rail[Double], dconA:Int, dconB:Int, temp:Rail[Double], N:Int, Ln:Int, Ylm:Rail[Double],maxL:Int):Int;

    @Native("c++","(#this)->getNL((#1)->raw)")
    public native def getNL(n_l:Rail[Int]):Int;
}

