package au.edu.anu.fft;

/**
 * Implements regular (non-fast) Discrete Fourier Transforms.
 */
public class DFT {
	public static def dft1D(input : Rail[Complex]!, forward : boolean) : Rail[Complex]! {
        val N = input.length;
        val output = Rail.make[Complex](N, (Int) => Complex.ZERO);
        for (var m:int=0; m<output.length; m++) {
            var sum : Complex = Complex.ZERO;        
            for (var k:int=0; k<input.length; k++) {
                output(m) = output(m) + input(k) * Complex.exp(Complex.I * 2.0 * Math.PI * ((m * k) / N)); 
            }
        }
        return output;
    }
}
