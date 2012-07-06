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
package edu.utk.cs.papi;

import x10.compiler.Native;
import x10.compiler.NativeRep;
import x10.compiler.NativeCPPCompilationUnit;
import x10.compiler.NativeCPPInclude;
import x10.util.IndexedMemoryChunk;

/**
 * C++ Wrapper for PAPI
 *
 * @author milthorpe
 */
@NativeCPPCompilationUnit("PAPI.cc")
@NativeRep("c++", "::edu::utk::cs::papi::PAPI *", "::edu::utk::cs::papi::PAPI", null)
public class PAPI {
    /**
     * Call to initialise PAPI.
     */
    @Native("c++", "(#this)->initialize()")
    public native def initialize():void;

    /**
     * Call to start PAPI FLOPS counters.
     */
    @Native("c++", "(#this)->startFlops()")
    public native def startFlops():void;

    /**
     * Call to stop currently running counters.
     */
    @Native("c++", "(#this)->stop()")
    public native def stop():void;

    /**
     * Get the value of a given counter.
     */
    @Native("c++", "(#this)->getCounter(#i)")
    public native def getCounter(i:Int):Long;   

}

