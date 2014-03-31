/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Josh Milthorpe 2011-2014.
 */
import x10.compiler.Inline;
import x10.io.File;
import x10.regionarray.Array;
import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;

/**
 * Benchmarks loop iteration using a variety 
 * of activity generation approaches
 * @author milthorpe 08/2011
 */
public class LocalityDistParallelLoop(size:Long, print:Boolean) {
    private static ITERS = 10;

    private static val work = new Rail[Long](Runtime.MAX_THREADS);

    public def this(size:Long, print:Boolean) {
        property(size, print);
    }

	public def testAll() {
        val a = DistArray.make[Long](Dist.makeBlock(Region.make(0..(size-1))));

        var start:Long = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                for ([p] in aLocal) async {
                    aLocal(p) = Runtime.workerId();
                    //work(Runtime.workerId())++;
                }
            }
        }
        var stop:Long = System.nanoTime();
        Console.OUT.printf("standard loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "standard.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                val regionIter = aLocal.region.iterator();
                while(regionIter.hasNext()) {
                    val p = regionIter.next();
                    async {
                        aLocal(p) = Runtime.workerId();
                        //work(Runtime.workerId())++;
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
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                val assign = (p:Long) => { aLocal(p) = Runtime.workerId(); /*work(Runtime.workerId())++;*/};
                LocalityDistParallelLoop.doForEach(aLocal.region, assign);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "bisect.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                val assign = (p:Long) => { aLocal(p) = Runtime.workerId(); /*work(Runtime.workerId())++;*/};
                finish LocalityDistParallelLoop.doForEach2(aLocal.region, assign);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("partially bisected loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "part_bisect.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                val nthreads = Runtime.NTHREADS;
                val chunkSize = size / nthreads;
                val remainder = size % nthreads;
                for (t in 0..(nthreads-1)) async {
                    val begin = (t < remainder) ? t*(chunkSize+1) : remainder + t*chunkSize;
                    val end = (t < remainder) ? (t+1)*(chunkSize+1) - 1 : remainder + (t+1)*chunkSize - 1;
                    for (p in begin..end) async {
                        aLocal(p) = Runtime.workerId();
                        //work(Runtime.workerId())++;
                    }
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("chunked loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "chunked.dat");
        printAndResetWork();

        start = System.nanoTime();
        for (i in 1..ITERS) {
            finish for (place in a.dist.places()) at(place) async {
                val aLocal = a.getLocalPortion() as Array[Long](1){rect};
                for ([p] in aLocal) {
                    aLocal(p) = Runtime.workerId();
                    //work(Runtime.workerId())++;
                }
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("sequential loop: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
        writeLocalSchedulingFile(a, "sequential.dat");
        printAndResetWork();
	}

    private def writeLocalSchedulingFile(a:DistArray[Long](1), fileName:String) {
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

    static class RegionBisection(r:Region) {
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
        public final def size() = r.size();
    }

    @Inline static final def doForEach(r:Region(1){rect==false}, closure:(Point) => void) {
        bisection(new RegionBisection(r), closure);
    }

    private static final def bisection(b:RegionBisection, closure:(Point) => void) {
        val size = b.size();
        if (size == 1) {
            val r = b.r;
            closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.max(i))));
        } else if (size == 2) {
            val r = b.r;
            async closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.min(i))));
            closure(Point.make(new Rail[Long](r.rank, (i:Long)=>r.max(i))));

        } else {
            val firstHalf = new RegionBisection(b.firstHalf);
            val secondHalf = new RegionBisection(b.secondHalf);
            async bisection(firstHalf, closure);
            bisection(secondHalf, closure);

        }
    }

    @Inline static final def doForEach(r:Region(1){rect}, closure:(Long) => void) {
        bisection1D(r.min(0), r.max(0), closure);
    }

    private static final def bisection1D(start:Long, end:Long, closure:(Long) => void) { 
        if (start == end) {
            closure(start);
        } else {
            val numElems = end-start+1;
            val halfway = start+numElems/2;
            async bisection1D(halfway, end, closure);
            bisection1D(start, halfway-1, closure);
        }
    }

    /** @return the depth of the binary tree spanning the threads */
    private static def getTreeDepth(numThreads:Int):Int {
        var p:Int = numThreads;
        var i:Int = 0n;
        while (p > 0n) { p = p>>1; i++; }
        return i;
    }

    @Inline final static def doForEach2(r:Region(1){rect}, closure:(Long) => void) {
        val treeDepth = getTreeDepth(Runtime.NTHREADS);
        bisection1D_limited(r.min(0), r.max(0), treeDepth, 0n, closure);
    }

    private static final def bisection1D_limited(start:Long, end:Long, treeDepth:Int, currentLevel:Int, closure:(Long) => void) {
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
            async bisection1D_limited(halfway, end, treeDepth, currentLevel+1n, closure);
            bisection1D_limited(start, halfway-1, treeDepth, currentLevel+1n, closure);
        }
    }

	public static def main(args:Rail[String]):void = {
        var size:Long = 8000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Long.parse(args(0));
            if (args.size > 1) {
                print=true;
            }
        }
		new LocalityDistParallelLoop(size, print).testAll();
	}

}
