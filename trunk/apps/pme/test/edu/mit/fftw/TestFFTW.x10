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
package edu.mit.fftw;

import x10.array.BaseArray;

public class TestFFTW {
    public static def main(args : Rail[String]!) {
        Console.OUT.println("1D");
        val N : Int = 50;
        val input = Rail.make[Complex](N, (i : Int) => 1.0 / (i+1.0) * (1.0 + Complex.I));
        for (var i : int = 0; i < input.length; i++) {
            Console.OUT.println("input(" + i + ") = " + input(i));
        }
        val output = Rail.make[Complex](N);
        val plan : FFTW.FFTWPlan = FFTW.fftwPlan1d(N, input, output, false);
        FFTW.fftwExecute(plan);
        FFTW.fftwDestroyPlan(plan);
        for (var i : int = 0; i < output.length; i++) {
            Console.OUT.println("output(" + i + ") =" + output(i));
        }
        
        val roundtrip = Rail.make[Complex](N);
        val plan2 : FFTW.FFTWPlan = FFTW.fftwPlan1d(N, output, roundtrip, true);
        FFTW.fftwExecute(plan2);
        FFTW.fftwDestroyPlan(plan2);
        for (var i : int = 0; i < roundtrip.length; i++) {
            Console.OUT.println("roundtrip(" + i + ") =" + roundtrip(i) / N);
        }

        Console.OUT.println("3D");
        val K = 4;
        val input3d = Array.make[Complex]([0..K-1, 0..K-1, 0..K-1], ((i,j,k) : Point) => 1.0 / (i+j+k+1.0) * (1.0 + Complex.I));
        val output3d =  Array.make[Complex](input3d.region);
        for (m in input3d) {
            Console.OUT.println("intput" + m + " = " + input3d(m));
        }
        val plan3 : FFTW.FFTWPlan = FFTW.fftwPlan3d(K, K, K, input3d as BaseArray[Complex], output3d as BaseArray[Complex], false);
        FFTW.fftwExecute(plan3);
        for (m in output3d) {
            Console.OUT.println("output" + m + " = " + output3d(m));
        }
        val roundtrip3d =  Array.make[Complex](input3d.region);
        val plan4 : FFTW.FFTWPlan = FFTW.fftwPlan3d(K, K, K, output3d as BaseArray[Complex], roundtrip3d as BaseArray[Complex], true);
        FFTW.fftwExecute(plan4);
        for (m in roundtrip3d) {
            Console.OUT.println("roundtrip" + m + " =" + roundtrip3d(m) / (K*K*K));
        }
    }
}
