package au.edu.anu.fft;

/**
 * Implements regular (non-fast) Discrete Fourier Transforms.
 */
public class DFT {
	public static def dft1D(input : Rail[Complex]!, forward : boolean) : Rail[Complex]! {
        val sign = forward ? -1.0 : 1.0;
        val n = input.length;
        val N = n as Double;
        val output = Rail.make[Complex](n, (Int) => Complex.ZERO);
        for (var m:int=0; m<n; m++) {
            val mDivN = (m as Double)* sign * 2.0 * Math.PI / N ;
            var mKDivN : Double = 0.0;
            var sum : Complex = Complex.ZERO;        
            for (var k:int=0; k<n; k++) {
                sum = sum + input(k) * Complex.exp(Complex.I * mKDivN); 
                mKDivN += mDivN;
            }
            output(m) = sum;
        }
        return output;
    }
}
