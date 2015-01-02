/*
 *  This file is part of the X10 project (http://x10-lang.org).
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2013.
 */

import x10.compiler.Pragma;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;

/**
 * Standard "ping/pong" benchmark for X10 communication performance.
 * This benchmark operates in two phases: the first uses remote array
 * copies and the second uses active messages/remote activities.
 * @author milthorpe 02/13
 */
public class Pong {
    static ITERS=100;

	public static def main(args:Rail[String]) {
        Console.OUT.println("Ping-pong using Remote Array Copy...");
        val dsrc = DistArray.make[Rail[Double]](Dist.makeUnique());
        val ddst = DistArray.make[Rail[Double]](Dist.makeUnique());
        val dremoteSrc = DistArray.make[Rail[GlobalRail[Double]]](Dist.makeUnique(), (Point)=> new Rail[GlobalRail[Double]](Place.numPlaces()));
        val dremoteDst = DistArray.make[Rail[GlobalRail[Double]]](Dist.makeUnique(), (Point)=> new Rail[GlobalRail[Double]](Place.numPlaces()));
        for (var s:Long=8; s<=8388608; s *= 2) {
            val size = s;
            // set up arrays at each place
            finish ateach(place in Dist.makeUnique()) {
                dsrc(here.id) = new Rail[Double](size, here.id as Double);
                ddst(here.id) = new Rail[Double](size);
            }
            // set up remote handles from each place to arrays at all other places
            finish ateach(place in Dist.makeUnique()) {
                for (place2 in Place.places()) {
                    dremoteSrc(here.id)(place2.id) = at(place2) new GlobalRail[Double](dsrc(here.id));
                    dremoteDst(here.id)(place2.id) = at(place2) new GlobalRail[Double](ddst(here.id));
                }
            }
            // ping-pong using remote handles
            // TODO ping-pong between all places, not just first two
            val nextId = Place.places().next(here).id;
            val localSrc = dsrc(here.id);
            val remoteSrc = dremoteSrc(here.id)(nextId);
            val localDst = ddst(here.id);
            val remoteDst = dremoteDst(here.id)(nextId);
            val start = System.nanoTime();
            for (i in 1..ITERS) {
                finish Rail.asyncCopy(localSrc, 0L, remoteDst, 0L, localSrc.size);
                finish Rail.asyncCopy(remoteSrc, 0L, localDst, 0L, localDst.size);
            }
            val stop = System.nanoTime();
            for (i in 0..(localDst.size-1)) {
                if (localDst(i) != nextId as Double) {
                    throw new Exception("incorrect value for element " + i + ": " + localDst(i));
                }
            }

            at(Place.places().next(here)) {
                val nextDst = ddst(here.id);
                val prevId = Place.places().prev(here).id;
                for (i in 0..(nextDst.size-1)) {
                    if (nextDst(i) != prevId as Double) {
                        throw new Exception("at " + here + " incorrect value for element " + i + ": " + nextDst(i));
                    }
                }
            }

            val time = ((stop-start) as Double) / 1e6 / ITERS;
            Console.OUT.printf("%d bytes took: %.3g ms (%.3g MB/sec)\n", size, time, 2.0*size / time / 1000);
        }

        Console.OUT.println("Ping-pong using active messages...");
        for (var size:Long=8; size<=8388608; size *= 2) {
            val a = new Rail[Double](size / 8);
            val start = System.nanoTime();
            for (i in 1..ITERS) {
                val b = at(Place.places().next(here)) {
                    a(0) = 1.0;
                    a
                };
            }
            val stop = System.nanoTime();
            val time = ((stop-start) as Double) / 1e6 / ITERS;
            Console.OUT.printf("%d bytes took: %.3g ms (%.3g MB/sec)\n", size, time, 2.0*size / time / 1000);
        }
	}
}
