/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
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

