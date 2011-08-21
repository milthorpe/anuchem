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

/**
 * Benchmarks loop iteration using a variety 
 * of activity generation approaches
 * @author milthorpe 08/2011
 */
public class BenchmarkParallelLoop(size : Int, print:Boolean) extends x10Test {
    private static ITERS = 1000;

    public def this(size : Int, print:Boolean) {
        property(size, print);
    }

	public def run(): Boolean = {
        val a = new Array[Int](size);

        var start:Long = System.nanoTime();
        for (i in 1..ITERS) {
            finish for ([p] in a) async {
                a(p) = p;
            }
        }
        var stop:Long = System.nanoTime();
        Console.OUT.printf("standard loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish {
                val regionIter = a.region.iterator();
                while(regionIter.hasNext()) {
                    val p = regionIter.next();
                    async {
                        a(p) = p(0);
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("iterator loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            val assign = (p:Int) => { a(p) = p;};
            finish BenchmarkParallelLoop.doForEach(a.region, assign);
        }
        stop = System.nanoTime();
        Console.OUT.printf("bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            val nthreads = Runtime.NTHREADS;
            val chunkSize = size / nthreads;
            val remainder = size % nthreads;
            finish for (t in 0..(nthreads-1)) async {
                val begin = (t < remainder) ? t*(chunkSize+1) : remainder + t*chunkSize;
                val end = (t < remainder) ? (t+1)*(chunkSize+1) - 1 : remainder + (t+1)*chunkSize - 1;
                for (p in begin..end) async {
                    a(p) = p;
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("chunked loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            for ([p] in a) {
                a(p) = p;
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("sequential loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        return true;
	}

    @Inline static def remainder(iter:Iterator[Point(1)], closure:(Point(1)) => void) {
        if (iter.hasNext()) {
            val p = iter.next();
            async remainder(iter, closure);
            closure(p);
        }
    }

    static struct RegionBisection(r:Region) {
        public val firstHalf:Region(r.rank);
        public val secondHalf:Region(r.rank);

        public def this(r:Region) {
            property(r);
            var firstHalf:Region(r.rank) = null;
            var secondHalf:Region(r.rank) = null;
            for (dim in 0..(r.rank-1)) {
                val min = r.min(dim);
                val max = r.max(dim);
                if (max > min) {
                    // bisect on this dimension
                    val halfway = (max-min+1) / 2;
                    val firstMax = new Array[Int](r.rank, (i:Int)=> i==dim ? (halfway-1) : r.max(i));
                    val firstCut = Region.makeRectangular(new Array[Int](r.rank, (i:Int)=>r.min(i)), firstMax);
                    firstHalf = r && firstCut;
                    val secondMin = new Array[Int](r.rank, (i:Int)=> i==dim ? halfway : r.min(i));
                    val secondCut = Region.makeRectangular(secondMin, new Array[Int](r.rank, (i:Int)=>r.max(i)));
                    secondHalf = r && secondCut;
                }
            }
            if (firstHalf != null) {
                this.firstHalf = firstHalf;
                this.secondHalf = secondHalf;
            } else {
                throw new IllegalArgumentException("unable to bisect region " + r);
            }
        }
        public def size() = r.size();
    }

    @Inline static final def doForEach(r:Region(1){rect==false}, closure:(Point) => void) {
        bisection(new RegionBisection(r), closure);
    }

    private static final def bisection(b:RegionBisection, closure:(Point) => void) {
        val size = b.size();
        if (size == 1) {
            val r = b.r;
            closure(Point.make(new Array[Int](r.rank, (i:Int)=>r.max(i))));
        } else if (size == 2) {
            val r = b.r;
            async closure(Point.make(new Array[Int](r.rank, (i:Int)=>r.max(i))));
            closure(Point.make(new Array[Int](r.rank, (i:Int)=>r.min(i))));
        } else {
            val firstHalf = new RegionBisection(b.firstHalf);
            val secondHalf = new RegionBisection(b.secondHalf);
            async bisection(secondHalf, closure);
            bisection(firstHalf, closure);
        }
    }

    @Inline final static def doForEach(r:Region(1){rect}, closure:(Int) => void) {
        bisection1D(r.min(0), r.max(0), closure);
    }

    private static final def bisection1D(start:Int, end:Int, closure:(Int) => void) { 
        if (start == end) {
            closure(start);
        } else {
            val numElems = end-start+1;
            val halfway = start+numElems/2;
            async bisection1D(halfway, end, closure);
            bisection1D(start, halfway-1, closure);
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
		new BenchmarkParallelLoop(size, print).execute();
	}

}
