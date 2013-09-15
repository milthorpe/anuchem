/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Josh Milthorpe 2011-2012.
 */
import x10.compiler.Inline;
import x10.regionarray.Array;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;

/**
 * Benchmarks dist array iteration using a variety 
 * of activity generation approaches
 * @author milthorpe 08/2011
 */
public class BenchmarkDistParallelLoop(size:Long, print:Boolean) {
    private static ITERS = 100;

    public def this(size:Long, print:Boolean) {
        property(size, print);
    }

	public def testAll() {
        val a = DistArray.make[Long](Dist.makeBlock(Region.make(0..(size-1))));

        var start:Long = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1);
                for ([p] in aLocal) async {
                    aLocal(p) = p;
                }
            }
        }
        var stop:Long = System.nanoTime();
        Console.OUT.printf("standard loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1);
                val regionIter = aLocal.region.iterator();
                while(regionIter.hasNext()) {
                    val p = regionIter.next();
                    async {
                        aLocal(p) = p(0);
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("iterator loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1);
                val assign = (p:Long) => { aLocal(p) = p;};
                BenchmarkDistParallelLoop.doForEach(aLocal.region, assign);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                val assign = (p:Long) => { aLocal(p) = p;};
                finish BenchmarkDistParallelLoop.doForEach2(aLocal.region, assign);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("partially bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1);
                val localSize = aLocal.size;
                val nthreads = Runtime.NTHREADS;
                val chunkSize = localSize / nthreads;
                val remainder = localSize % nthreads;
                for (t in 0..(nthreads-1)) async {
                    val begin = aLocal.region.min(0) + ((t < remainder) ? t*(chunkSize+1) : remainder + t*chunkSize);
                    val end = aLocal.region.min(0) + ((t < remainder) ? (t+1)*(chunkSize+1) - 1 : remainder + (t+1)*chunkSize - 1);
                    for (p in begin..end) async {
                        aLocal(p) = p;
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("chunked loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1);
                for ([p] in aLocal) {
                    aLocal(p) = p;
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("sequential loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
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
                    val firstMax = new Rail[Long](r.rank, (i:Long)=> i==dim ? (halfway-1) : r.max(i));
                    val firstCut = Region.makeRectangular(new Rail[Long](r.rank, (i:Long)=>r.min(i)), firstMax);
                    firstHalf = r && firstCut;
                    val secondMin = new Rail[Long](r.rank, (i:Long)=> i==dim ? halfway : r.min(i));
                    val secondCut = Region.makeRectangular(secondMin, new Rail[Long](r.rank, (i:Long)=>r.max(i)));
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

    @Inline static def doForEach(r:Region(1){rect==false}, closure:(Point) => void) {
        bisection(new RegionBisection(r), closure);
    }

    private static def bisection(b:RegionBisection, closure:(Point) => void) {
        val size = b.size();
        if (size == 1) {
            val r = b.r;
            closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.max(i))));
        } else if (size == 2) {
            val r = b.r;
            async closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.max(i))));
            closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.min(i))));
        } else {
            val firstHalf = new RegionBisection(b.firstHalf);
            val secondHalf = new RegionBisection(b.secondHalf);
            async bisection(secondHalf, closure);
            bisection(firstHalf, closure);
        }
    }

    @Inline static def doForEach(r:Region(1){rect}, closure:(Long) => void) {
        bisection1D(r.min(0), r.max(0), closure);
    }

    private static def bisection1D(start:Long, end:Long, closure:(Long) => void) { 
        val numElems = end-start+1;
        if (numElems == 1) {
            closure(start);
        } else if (numElems == 2) {
            async closure(end);
            closure(start);
        } else {
            val halfway = start+numElems/2;
            async bisection1D(halfway, end, closure);
            bisection1D(start, halfway-1, closure);
        }
    }

    /** @return the depth of the binary tree spanning the threads */
    private static def getTreeDepth(numThreads:Long):Long {
        var p:Long = numThreads;
        var i:Long = 0;
        while (p > 0) { p = p>>1; i++; }
        return i;
    }

    @Inline final static def doForEach2(r:Region(1){rect}, closure:(Long) => void) {
        val treeDepth = getTreeDepth(Runtime.NTHREADS);
        bisection1D_limited(r.min(0), r.max(0), treeDepth, 0, closure);
    }

    private static final def bisection1D_limited(start:Long, end:Long, treeDepth:Long, currentLevel:Long, closure:(Long) => void) {
        if (currentLevel == treeDepth) {
            // don't divide any further
            for (p in start..end) async {
                closure(p);

            }
        } else if (start == end) {
            closure(start);
        } else {
            val numElems = end-start+1;
            val halfway = start+numElems/2;
            async bisection1D_limited(halfway, end, treeDepth, currentLevel+1, closure);
            bisection1D_limited(start, halfway-1, treeDepth, currentLevel+1, closure);
        }
    }

	public static def main(var args:Rail[String]): void = {
        var size:Long = 8000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Long.parse(args(0));
            if (args.size > 1) {
                print=true;
            }
        }
		new BenchmarkDistParallelLoop(size, print).testAll();
	}
}
