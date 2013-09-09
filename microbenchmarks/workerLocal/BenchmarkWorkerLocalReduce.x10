/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright Australian National University 2013.
 */
import x10.compiler.Inline;
import x10.util.WorkerLocalHandle;

/**
 * Benchmarks reducing worker local data
 * @author milthorpe 06/2013
 */
public class BenchmarkWorkerLocalReduce(N:Long, print:Boolean) {
    private static ITERS = 1000;

    public def this(N:Long, print:Boolean) {
        property(N, print);
    }

	public def testAll() {
        var start:Long;
        var stop:Long;
/*
        val result_worker = new WorkerLocalHandle[Rail[Double]]( 
            () => new Rail[Double](N)
        );
        finish for(i in 0..(N*Runtime.NTHREADS-1)) async {
            val partialResult = result_worker();
            partialResult(i%N) += i;
        }
        var result:Rail[Double] = null;

        start = System.nanoTime();
        for (i in 1..ITERS) {
            result = result_worker.reduceLocal( 
                (a:Rail[Double],b:Rail[Double]) => a.map(a, b, 
                    (x:Double,y:Double)=>(x+y)) as Rail[Double]
            );
        }
        stop = System.nanoTime();
        val x = result.reduce((a:Double,b:Double)=>(a+b), 0.0);
        Console.OUT.println(x);
        Console.OUT.printf("WLH.reduceLocal rail: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        // compare with cost of reducing Runtime.NTHREADS partial result arrays in sequence
        val partial = new Rail[Rail[Double]](Runtime.NTHREADS, (t:Int) => new Rail[Double](N, (i:Int)=>(i+t*N) as Double));

        var result2:Rail[Double] = null;
        start = System.nanoTime();
        for (i in 1..ITERS) {
            result2 = new Rail[Double](N);
            for (j in 0..(Runtime.NTHREADS-1)) {
                result2.map(result2, partial(j),
                        (x:Double,y:Double)=>(x+y)
                );
            }
        }
        stop = System.nanoTime();
        val y = result2.reduce((a:Double,b:Double)=>(a+b), 0.0);
        Console.OUT.println(y);
        Console.OUT.printf("sequential reduce rail: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
*/
        var sum:Rail[Double] = null;
        start = System.nanoTime();
        for (i in 1..ITERS) {
            val sum_worker = new WorkerLocalHandle[Rail[Double]]( 
                () => new Rail[Double](1L)
            );
            finish for (j in 1..N) async {
                sum_worker()(0)++;
            }
            sum = sum_worker.reduceLocal( 
                (a:Rail[Double],b:Rail[Double]) => {
                    a(0) += b(0);
                    a
                }
            );
        }
        stop = System.nanoTime();
        Console.OUT.println(sum(0));
        Console.OUT.printf("WLH reduce double: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        var sum2:Double = 0.0;
        start = System.nanoTime();
        for (i in 1..ITERS) {
            sum2 = finish(Reducible.SumReducer[Double]()) {
                for (j in 1..N) async {
                    offer 1.0;
                }
            };
        }
        stop = System.nanoTime();
        Console.OUT.println(sum2);
        Console.OUT.printf("collecting finish double: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
	}

	public static def main(args:Rail[String]): void = {
        var size:Long = 100000;
        var print:Boolean = false;
        if (args.size > 0) {
            size = Long.parse(args(0));
            if (args.size > 1) {
                print=true;
            }
        }
		new BenchmarkWorkerLocalReduce(size, print).testAll();
	}
}
