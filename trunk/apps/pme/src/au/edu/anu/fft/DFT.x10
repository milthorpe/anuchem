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
                sum = sum + input(k) * (Complex.I * mKDivN).exp(); 
                mKDivN += mDivN;
            }
            output(m) = sum;
        }
        return output;
    }

    public static def dft1D(input : Array[Complex](3)!{self.rect,self.zeroBased}, forward : boolean) : Array[Complex](3)! {
        val sign = forward ? -1.0 : 1.0;
        val r = input.region;
        val n1 = r.max(0) + 1;
        val n2 = r.max(1) + 1;
        val n3 = r.max(2) + 1;
        val N1 = n1 as Double;
        val N2 = n2 as Double;
        val N3 = n3 as Double;
        val output = Array.make[Complex](input.region);
        for (var m1:int=0; m1<n1; m1++) {
            for (var m2:int=0; m2<n2; m2++) {
                for (var m3:int=0; m3<n3; m3++) {
                    var sum : Complex = Complex.ZERO;        
                    for (var k1:int=0; k1<n1; k1++) {
                        for (var k2:int=0; k2<n2; k2++) {
                            for (var k3:int=0; k3<n3; k3++) {
                                val component = sign * 2.0 * Math.PI * Complex.I * ((m1 * k1 as Double) / N1 + (m2 * k2 as Double) / N2 + (m3 * k3 as Double) / N3);
                                sum = sum + input(k1, k2, k3) * component.exp();
                            }
                        }
                    }
                    output(m1, m2, m3) = sum;
                }
            }
        }
        return output;
    }
}
