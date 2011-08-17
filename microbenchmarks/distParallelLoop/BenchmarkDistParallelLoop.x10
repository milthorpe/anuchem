/*
 *  This file is part of the X10 project (http://x10-lang.org).
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2011.
 */
import harness.x10Test;

import x10.compiler.Inline;
import x10.io.File;

/**
 * Benchmarks loop iteration using a variety 
 * of activity generation approaches
 * @author milthorpe 06/2011
 */
public class BenchmarkDistParallelLoop(size : Int, print:Boolean) extends x10Test {
    private static ITERS = 10;

    private static val work = new Array[Int](Runtime.MAX_THREADS);

    public def this(size : Int, print:Boolean) {
        property(size, print);
    }

	public def run(): Boolean = {
        val a = DistArray.make[Int](Dist.makeBlock(0..(size-1)));

        var start:Long = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) async at(place) {
                val aLocal = a.getLocalPortion() as Rail[Int];
                for (p in aLocal) async {
                    aLocal(p) = Runtime.workerId();
                    work(Runtime.workerId())++;
                }
            }
        }
        var stop:Long = System.nanoTime();
        Console.OUT.printf("standard loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "standard.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) async at(place) {
                val aLocal = a.getLocalPortion() as Rail[Int];
                val regionIter = aLocal.region.iterator();
                while(regionIter.hasNext()) {
                    val p = regionIter.next();
                    async {
                        aLocal(p) = Runtime.workerId();
                        work(Runtime.workerId())++;
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("iterator loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "iterator.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) async at(place) {
                val aLocal = a.getLocalPortion();
                val regionIter = aLocal.region;
                val assign = (p:Int) => { aLocal(p) = Runtime.workerId(); work(Runtime.workerId())++;};
                BenchmarkDistParallelLoop.doBisectingLoop(regionIter.min(0), regionIter.max(0), assign);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "bisect.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) async at(place) {
                val aLocal = a.getLocalPortion();
                val nthreads = Runtime.NTHREADS;
                val chunkSize = size / nthreads;
                val remainder = size % nthreads;
                for (t in 0..(nthreads-1)) {
                    val begin = (t < remainder) ? t*(chunkSize+1) : remainder + t*chunkSize;
                    val end = (t < remainder) ? (t+1)*(chunkSize+1) - 1 : remainder + (t+1)*chunkSize - 1;
                    async for (p in begin..end) {
                        aLocal(p) = Runtime.workerId();
                        work(Runtime.workerId())++;
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("chunked loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "chunked.dat");
        printAndResetWork();

        return true;
	}

    private def writeLocalSchedulingFile(a:DistArray[Int](1), fileName:String) {
        val outputFile = new File(fileName);
        val out = outputFile.printer();
        val aLoc = a.getLocalPortion();
        for (p in aLoc) {
            out.println(p(0) + " " + aLoc(p));
        }
        out.flush();
    }

    def printAndResetWork() {
        for (w in work) {
            if (work(w) > 0) {
                Console.OUT.println("worker " + w + " = " + work(w));
            }
            work(w) = 0;
        }
    }

    @Inline static def remainder(iter:Iterator[Point(1)], closure:(Point(1)) => void) {
        if (iter.hasNext()) {
            val p = iter.next();
            async remainder(iter, closure);
            closure(p);
        }
    }

    @Inline static def doBisectingLoop(start:Int, end:Int, closure:(Int) => void) { 
        val numElems = end-start+1;
        if (numElems == 1) {
            closure(start);
        } else if (numElems == 2) {
            async closure(end);
            closure(start);
        } else {
            val halfway = start+numElems/2;
            async doBisectingLoop(halfway, end, closure);
            doBisectingLoop(start, halfway-1, closure);
        }
    }

	public static def main(var args: Array[String](1)): void = {
        var size:Int = 8000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Int.parse(args(0));
            if (args.size > 1) {
                print=true;
            }
        }
		new BenchmarkDistParallelLoop(size, print).execute();
	}

}
