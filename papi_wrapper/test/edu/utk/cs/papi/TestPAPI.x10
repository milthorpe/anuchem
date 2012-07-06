/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package edu.utk.cs.papi;

/**
 * Tests wrapper for PAPI toolkit
 * @see http://icl.cs.utk.edu/papi
 * @author milthorpe
 */
public class TestPAPI {
    public def sizeOfCentralCluster() : Double = 80.0;

    public static def main(args : Array[String](1)) {
        val papi = new PAPI();
        papi.initialize();
        papi.startFlops();
        papi.stop();
        Console.OUT.printf("Total cycles: %d total FP ops: %d", papi.getCounter(0), papi.getCounter(1));
    }
}

