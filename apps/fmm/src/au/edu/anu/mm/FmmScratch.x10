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
package au.edu.anu.mm;

/**
 * This class provides worker-local scratch space for use in the 
 * FMM operators, to avoid overhead of array allocation and GC.
 */
public class FmmScratch {
    var exp:MultipoleExpansion;
    var array:Rail[Complex];
    var plm:AssociatedLegendrePolynomial;

    public def this(numTerms:Int) {
        exp = new MultipoleExpansion(numTerms);
        array = new Array[Complex](numTerms+1);
        plm = new AssociatedLegendrePolynomial(numTerms);
    }

    private static val store:Rail[FmmScratch] = new Rail[FmmScratch](Runtime.MAX_THREADS);

    public static getWorkerLocal():FmmScratch {
        return store(Runtime.workerId());
    }

    public static def setWorkerLocal(scratch:FmmScratch):void {
        store(Runtime.workerId()) = scratch;
    }

    public static def init(func:()=>FmmScratch):void {
        for (i in 0..(store.size-1)) {
            store(i) = func();
        }
    }
}

