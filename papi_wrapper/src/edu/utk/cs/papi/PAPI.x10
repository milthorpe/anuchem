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
    static COUNTER_TOTAL_INS = 0;

    // floating point ops
    static COUNTER_FP_INS = 1;
    static COUNTER_FP_OPS = 2;

    // memory ops    
    static COUNTER_LD_INS = 1;
    static COUNTER_SR_INS = 2;
    static COUNTER_L1_LDM = 3;
    static COUNTER_L1_STM = 4;

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
     * Creates eventset for FLOPS counters.
     * Assumes counters are not currently running.
     */
    @Native("c++", "(#this)->countFlops()")
    public native def countFlops():void;

    /**
     * Print current values of FLOPS counters
     */
    @Native("c++", "(#this)->printFlops()")
    public native def printFlops():void;

    /**
     * Creates eventset for load/store counters.
     * Assumes counters are not currently running.
     */
    @Native("c++", "(#this)->countMemoryOps()")
    public native def countMemoryOps():void;

    /**
     * Print current values of load/store counters
     */
    @Native("c++", "(#this)->printMemoryOps()")
    public native def printMemoryOps():void;

    /**
     * Start or resume counters for current event set
     */
    @Native("c++", "(#this)->start()")
    public native def start():void;

    /**
     * Stop and accumulate currently running counters
     */
    @Native("c++", "(#this)->stop()")
    public native def stop():void;

    /**
     * Reset currently running counters
     */
    @Native("c++", "(#this)->reset()")
    public native def reset():void;

    /**
     * Get the value of a given counter
     */
    @Native("c++", "(#this)->getCounter(#i)")
    public native def getCounter(i:Int):Long;

    /**
     * Set the value of a given counter to zero
     */
    @Native("c++", "(#this)->resetCounter(#i)")
    public native def resetCounter(i:Int):void;   

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

