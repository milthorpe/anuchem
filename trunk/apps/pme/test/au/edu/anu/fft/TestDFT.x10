package au.edu.anu.fft;


/**
 * Tests regular (non-fast) Discrete Fourier Transforms.
 */
public class TestDFT {
    public static def main(args : Rail[String]!) {
        val N = 10;

        val onePlusI = Complex.ONE + Complex.I;
        val input = ValRail.make[Complex](N, (i : Int) => onePlusI * (1.0 / (i + 1)));
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
}
