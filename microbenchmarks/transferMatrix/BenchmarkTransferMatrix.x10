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

import x10.compiler.Pragma;
import x10.matrix.block.Grid;
import x10.matrix.DenseMatrix;
import x10.matrix.dist.DistDenseMatrix;

/**
 * Benchmarks retrieving the local portion of a DistDenseMatrix from another place
 * @author milthorpe 10/2013
 */
public class BenchmarkTransferMatrix(M:Long, N:Long) {
    private static ITERS = 100;

    public def this(M:Long, N:Long) {
        property(M, N);
    }

	public def test() {
        var start:Long;
        var stop:Long;

        val grid = Grid.make(M, N);
        val ddm = DistDenseMatrix.make(grid);

        var nextBlock:DenseMatrix = null;
        start = System.nanoTime();
        for (i in 1..ITERS) {
            nextBlock = at(here.next()) {ddm.local()};
        }
        stop = System.nanoTime();
        Console.OUT.println("nextBlock size = " + (nextBlock.M*nextBlock.N) + " ("+ nextBlock.M + ", " + nextBlock.N+")");
        Console.OUT.printf("transfer with at expression: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);

        // copy from next place to here using Rail.asyncCopy
        nextBlock = new DenseMatrix(ddm.local().M, ddm.local().N);
        val hereRef = new GlobalRail(ddm.local().d);
        start = System.nanoTime();
        for (i in 1..ITERS) {
            //nextBlock = at(here.next()) {ddm.local()};
            @Pragma(Pragma.FINISH_ASYNC_AND_BACK) finish at(here.next()) {
                val dataHere = ddm.local().d;
                Rail.asyncCopy(dataHere, 0, hereRef, 0, dataHere.size);
            }
        }
        stop = System.nanoTime();
        Console.OUT.printf("transfer with Rail.asyncCopy: %g ms\n", ((stop-start) as Double) / 1e6 / ITERS);
	}

	public static def main(args:Rail[String]): void {
        var m:Long = 14450;
        var n:Long = 700;
        var print:Boolean = false;
        if (args.size > 0) {
            m = Long.parse(args(0));
            if (args.size > 1) {
                n = Long.parse(args(1));
            }
        }
		new BenchmarkTransferMatrix(m, n).test();
	}
}
