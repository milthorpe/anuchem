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

import x10.util.IndexedMemoryChunk;
import x10.util.RemoteIndexedMemoryChunk;

/**
 * Standard "ping-pong" benchmark for X10 communication performance.
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
        val dremoteSrc = DistArray.make[Rail[RemoteIndexedMemoryChunk[Double]]](Dist.makeUnique());
        val dremoteDst = DistArray.make[Rail[RemoteIndexedMemoryChunk[Double]]](Dist.makeUnique());
        for (var s:Int=8; s<=8388608; s *= 2) {
            val size = s;
            // set up arrays at each place
            finish ateach(place in Dist.makeUnique()) {
                dsrc(here.id) = new Rail[Double](size, here.id as Double);
                ddst(here.id) = new Rail[Double](size);
            }
            // set up remote handles from each place to arrays at all other places
            finish ateach(place in Dist.makeUnique()) {
                for (place2 in PlaceGroup.WORLD) {
                    dremoteSrc(here.id)(place2.id) = at(place2) RemoteIndexedMemoryChunk.wrap[Double](dsrc(here.id).raw());
                    dremoteDst(here.id)(place2.id) = at(place2) RemoteIndexedMemoryChunk.wrap[Double](ddst(here.id).raw());
                }
            }
            // ping-pong using remote handles
            // TODO ping-pong between all places, not just first two
            val nextId = here.next().id;
            val localSrc = dsrc(here.id).raw();
            val remoteSrc = dremoteSrc(here.id)(nextId);
            val localDst = ddst(here.id).raw();
            val remoteDst = dremoteDst(here.id)(nextId);
            val start = System.nanoTime();
            for (i in 1..ITERS) {
                finish IndexedMemoryChunk.asyncCopy(localSrc, 0, remoteDst, 0, localSrc.length());
                finish IndexedMemoryChunk.asyncCopy(remoteSrc, 0, localDst, 0, localDst.length());
            }
            val stop = System.nanoTime();
            for (i in 0..(localDst.length()-1)) {
                if (localDst(i) != nextId as Double) {
                    throw new Exception("incorrect value for element " + i + ": " + localDst(i));
                }
            }

            at(here.next()) {
                val prevId = here.prev().id;
                for (i in 0..(localDst.length()-1)) {
                    if (localDst(i) != prevId as Double) {
                        throw new Exception("at next incorrect value for element " + i + ": " + localDst(i));
                    }
                }
            }

            val time = ((stop-start) as Double) / 1e6 / ITERS;
            Console.OUT.printf("%d bytes took: %.3g ms (%.3g MB/sec)\n", size, time, 2.0*size / time / 1000);
        }

        Console.OUT.println("Ping-pong using active messages...");
        for (var size:Int=8; size<=8388608; size *= 2) {
            val a = new Array[Double](size / 8);
            val start = System.nanoTime();
            for (i in 1..ITERS) {
                val b = at(here.next()) {
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
