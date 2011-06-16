/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.edu.anu.util;

import x10.util.concurrent.AtomicInteger;

/**
 * A simple shared counter.
 *
 * Usage: 
 *      val sc = new SharedCounter();
 *      finish for(place in Place.places) {
 *          async at(place) {
 *              val current = sc.increment();
 *          }
 *      }
 *      val final = sc.get();
 *
 * @author milthorpe
 */
public class SharedCounter {
    private var counter : GlobalRef[AtomicInteger];

    public def this() {
        val a = new AtomicInteger();
        counter = GlobalRef[AtomicInteger](a);
    }

    public def getAndIncrement() { 
        return at(counter) {(counter as GlobalRef[AtomicInteger]{self.home==here})().getAndIncrement()};
    }
 
    public def get() {
        return at(counter) {(counter as GlobalRef[AtomicInteger]{self.home==here})().get()};
    }

    public def set(v : Int) {
        at(counter) {(counter as GlobalRef[AtomicInteger]{self.home==here})().set(v);};
    }
}

