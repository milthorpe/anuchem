package x10x.matrix;

import x10.util.Random;
import harness.x10Test;

/**
 * Benchmarks XLA Matrix methods
 * @author milthorpe
 */
public class BenchmarkMatrix extends x10Test {
    val N : Int;
    val rand = new Random(47);

    public def this(N : Int) {
        super();
        this.N = N;
    }

    public def run(): boolean {
        Console.OUT.println("square matrix size " + N);
        val a = randomMatrix();
        val start = System.nanoTime();
        val x = a.sumOffDiagonal();
        val stop = System.nanoTime();
        Console.OUT.printf("sumOffDiagonal: %g ms\n", ((stop-start) as Double) / 1e6);

        val start2 = System.nanoTime();
        val a2 = a.mul(2.0);
        val stop2 = System.nanoTime();
        Console.OUT.printf("scale: %g ms\n", ((stop2-start2) as Double) / 1e6);

        val b = randomMatrix();
        val start3 = System.nanoTime();
        val c = a.mul(b);
        val stop3 = System.nanoTime();
        Console.OUT.printf("gemm: %g ms\n", ((stop3-start3) as Double) / 1e6);

        val start4 = System.nanoTime();
        val d = a.transpose();
        val stop4 = System.nanoTime();
        Console.OUT.printf("transpose: %g ms\n", ((stop4-start4) as Double) / 1e6);

        return true;
    }

    private def randomMatrix() {
        val a = new Matrix(N);
        for ([i,j] in a.region) {
            a(i,j) = rand.nextDouble();
        }
        return a;
    }

    public static def main(args : Array[String](1)) {
        var N : Int = 70;
        if (args.size > 0) {
            N = Int.parseInt(args(0));
        }
        new BenchmarkMatrix(N).execute();
    }
}
