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

package au.edu.anu.util;

import x10.array.DistArray;

public class SharedCounter {
   global val localCounters:DistArray[Double];

   public def this() {
       localCounters = DistArray.make[Double](Dist.makeUnique());
   }

   public global def increment() = localCounters(here.id) += 1;
 
   public global def increment(step:Double) = localCounters(here.id) += step;

   public def get() = localCounters.reduce((a:Double,b:Double)=>a+b, 0.0);
}

