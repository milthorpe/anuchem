/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
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
                sum = sum + input(k) * Math.exp(Complex.I * mKDivN); 
                mKDivN += mDivN;
            }
            output(m) = sum;
        }
        return output;
    }

    public static def dft3D(input : Array[Complex](3)!{self.rect,self.zeroBased}, forward : boolean) : Array[Complex](3)! {
        val sign = forward ? -1.0 : 1.0;
        val r = input.region;
        val n1 = r.max(0) + 1;
        val n2 = r.max(1) + 1;
        val n3 = r.max(2) + 1;
        val invN1 = 1.0 / (n1 as Double);
        val invN2 = 1.0 / (n2 as Double);
        val invN3 = 1.0 / (n3 as Double);

        val s = sign * 2.0 * Math.PI * Complex.I;

        val output = new Array[Complex](input.region);
        for (var m1:int=0; m1<n1; m1++) {
            for (var m2:int=0; m2<n2; m2++) {
                for (var m3:int=0; m3<n3; m3++) {
                    var sum : Complex = Complex.ZERO;        
                    for (var k1:int=0; k1<n1; k1++) {
                        for (var k2:int=0; k2<n2; k2++) {
                            for (var k3:int=0; k3<n3; k3++) {
                                val component = s * (m1 * k1 * invN1 + m2 * k2 * invN2 + m3 * k3 * invN3);
                                sum = sum + input(k1, k2, k3) * Math.exp(component);
                            }
                        }
                    }
                    output(m1, m2, m3) = sum;
                }
            }
        }
        return output as Array[Complex](3)!;
    }
}
