package x10x.matrix;

import x10.util.Random;
import harness.x10Test;

/**
 * Tests XLA Matrix class multiplication method
 * @author milthorpe
 */
public class TestMatrixMultiply extends x10Test {
    val N : Int;
    val rand = new Random(47);

    public def this(N : Int) {
        super();
        this.N = N;
    }

    public def run(): boolean {
        val a = randomMatrix();
        Console.OUT.println(a);
        val b = randomMatrix();
        Console.OUT.println(b);

        val c = a.mul(b);

        Console.OUT.println(c);

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
        var N : Int = 4;
        if (args.size > 0) {
            N = Int.parseInt(args(0));
        }
        new TestMatrixMultiply(N).execute();
    }
}
