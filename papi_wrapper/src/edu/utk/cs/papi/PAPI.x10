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
    // High-level API

    /**
     * Initialize PAPI
     */
    @Native("c++", "(#this)->initialize()")
    public native def initialize():void;

    /**
     * Destroy current eventset and shutdown PAPI
     */
    @Native("c++", "(#this)->shutDown()")
    public native def shutDown():void;

    /**
     * Add events for FLOPS counters and start count
     */
    @Native("c++", "(#this)->countFlops()")
    public native def countFlops():void;

    /**
     * Print FLOPS counters
     */
    @Native("c++", "(#this)->printFlops()")
    public native def printFlops():void;

    /**
     * Add events for load/store counters and start count
     */
    @Native("c++", "(#this)->countMemoryOps()")
    public native def countMemoryOps():void;

    /**
     * Print load/store counters
     */
    @Native("c++", "(#this)->printMemoryOps()")
    public native def printMemoryOps():void;

    /**
     * Start counters for current event set
     */
    @Native("c++", "(#this)->startCount()")
    public native def startCount():void;

    /**
     * Stop currently running counters
     */
    @Native("c++", "(#this)->stopCount()")
    public native def stopCount():void;

    /**
     * Resume counters for current event set
     */
    @Native("c++", "(#this)->resumeCount()")
    public native def resumeCount():void;

    /**
     * Reset currently running counters and leave running
     */
    @Native("c++", "(#this)->resetCount()")
    public native def resetCount():void;

    /**
     * Get the value of a given counter
     */
    @Native("c++", "(#this)->getCounter(#i)")
    public native def getCounter(i:Int):Long;   

    // Low-level API

    /**
     * Create an eventset for PAPI counters
     */
    @Native("c++", "(#this)->createEventSet()")
    public native def createEventSet():void;

    /**
     * Add an event to the current eventset
     */
    @Native("c++", "(#this)->addEvent(#eventCode)")
    public native def addEvent(eventCode:Int):void;

    /**
     * Destroy current eventset
     */
    @Native("c++", "(#this)->destroyEventSet()")
    public native def destroyEventSet():void;

}

