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

import x10.array.DistArray;

/**
 * A simple shared counter.
 *
 * Usage: 
 *      val sc = new SharedCounter();
 *      finish for(place in Place.places) {
 *          async at(place) {
 *              sc.increment();
 *          }
 *      }
 *      val total = sc.get();
 *      Console.OUT.println(total);
 *      Console.OUT.println("success? : " + (total == Place.places.length));
 *
 * @author V. Ganesh, with modifications from Igor
 */
public class SharedCounter {
   global val localCounters:DistArray[Double];

   public def this() {
       localCounters = DistArray.make[Double](Dist.makeUnique());
   }

   public global def increment() = localCounters(here.id) += 1;
 
   public global def increment(step:Double) = localCounters(here.id) += step;

   public def get() = localCounters.reduce((a:Double,b:Double)=>a+b, 0.0);
}

