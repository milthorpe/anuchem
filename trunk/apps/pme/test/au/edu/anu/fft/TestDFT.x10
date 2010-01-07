package au.edu.anu.fft;


/**
 * Tests regular (non-fast) Discrete Fourier Transforms.
 */
public class TestDFT {
    public static def main(args : Rail[String]!) {
        testDft1D();
        testDft3D();
    }

    public static def testDft1D() {
        val N = 10;

        val onePlusI = Complex.ONE + Complex.I;
        val input = Rail.make[Complex](N, (i : Int) => onePlusI * (1.0 / (i + 1))) as Rail[Complex]!;
        Console.OUT.println("input");
        for (var i:int=0; i<input.length; i++) {
            Console.OUT.println(input(i));
        }

        Console.OUT.println("\noutput");
        val output = DFT.dft1D(input, false);
        for (var i:int=0; i<output.length; i++) {
            Console.OUT.println(output(i));
        }

        Console.OUT.println("\nroundtrip");
        val roundtrip = DFT.dft1D(output, true);
        for (var i:int=0; i<roundtrip.length; i++) {
            Console.OUT.println(roundtrip(i) / N);
        }
    }

    public static def testDft3D() {
        val N = 3;

        val twoPlusI = 2.0 + Complex.I;
        val r = Region.make(0, N-1);
        val r3 = r * r * r as Region(3);
        val input = Array.make[Complex](r3, (p(i,j,k) : Point(3)) => twoPlusI * (1.0 / (i + j + k + 1)));
        Console.OUT.println("input");
        for (p in input) {
            Console.OUT.println(input(p));
        }

        Console.OUT.println("\noutput");
        val output = DFT.dft1D(input, false);
        for (p in output) {
            Console.OUT.println(output(p));
        }

        Console.OUT.println("\nroundtrip");
        val roundtrip = DFT.dft1D(output, true);
        for (p in roundtrip) {
            Console.OUT.println(roundtrip(p) / (N*N*N));
        }
    }
}
