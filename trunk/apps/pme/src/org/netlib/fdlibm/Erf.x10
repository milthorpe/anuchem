package org.netlib.fdlibm;

import x10.compiler.Inline;
import x10.compiler.Native;
import x10.compiler.NativeCPPCompilationUnit;
import x10.compiler.NativeCPPInclude;

@NativeCPPInclude("erf.h")
@NativeCPPCompilationUnit("erf.cc")
public class Erf {
    @Native("c++", "fdlibm_erf(#1)")
    public native static @Inline def erf(x : Double) : Double;

    @Native("c++", "fdlibm_erfc(#1)")
    public native static @Inline def erfc(x : Double) : Double;
}

