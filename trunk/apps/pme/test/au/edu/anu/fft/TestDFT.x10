package au.edu.anu.fft;


/**
 * Tests regular (non-fast) Discrete Fourier Transforms.
 */
public class TestDFT {
    public static def main(args : Rail[String]!) {
        val onePlusI = Complex.ONE + Complex.I;
        val input = ValRail.make[Complex](50, (i : Int) => onePlusI * (1.0 / (i + 1)));
        val output = DFT.dft1D(input, true);
        for (var i:int=0; i<output.length; i++) {
            Console.OUT.println(output(i));
        }
    }
}
