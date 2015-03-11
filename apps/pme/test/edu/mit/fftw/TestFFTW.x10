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

import x10.regionarray.Dist;
import x10.regionarray.DistArray;
import x10.regionarray.Region;

public class TestFFTW {
    public static def main(args:Rail[String]) {
        Console.OUT.println("1D");
        val N:Int = 50n;
        val input = new Rail[Complex](N, (i:Long) => 1.0 / (i+1.0) * (1.0 + Complex.I));
        for (i in 0..(input.size-1)) {
            Console.OUT.println("input(" + i + ") = " + input(i));
        }
        val output = new Rail[Complex](N);
        val plan = FFTW.fftwPlan1d(N, input, output, false);
        FFTW.fftwExecute(plan);
        FFTW.fftwDestroyPlan(plan);
        for (i in 0..(output.size-1)) {
            Console.OUT.println("output(" + i + ") =" + output(i));
        }
        
        val roundtrip = new Rail[Complex](N);
        val plan2 = FFTW.fftwPlan1d(N, output, roundtrip, true);
        FFTW.fftwExecute(plan2);
        FFTW.fftwDestroyPlan(plan2);
        for (i in 0..(roundtrip.size-1)) {
            Console.OUT.println("roundtrip(" + i + ") =" + roundtrip(i) / N);
        }

        Console.OUT.println("3D");
        val K = 4n;
        val r = Region.makeRectangular(0..(K-1), 0..(K-1), 0..(K-1));
        val d = Dist.makeBlockBlock(r, 0, 1);
        val input3d = DistArray.make[Complex](d, ([i,j,k] : Point) => 1.0 / (i+j+k+1.0) * (1.0 + Complex.I));
        val output3d = DistArray.make[Complex](d);
        for (m in input3d) {
            Console.OUT.println("intput" + m + " = " + input3d(m));
        }
        val plan3 = FFTW.fftwPlan3d(K, K, K, input3d, output3d, false);
        FFTW.fftwExecute(plan3);
        for (m in output3d) {
            Console.OUT.println("output" + m + " = " + output3d(m));
        }
        val roundtrip3d =  DistArray.make[Complex](d);
        val plan4 = FFTW.fftwPlan3d(K, K, K, output3d, roundtrip3d, true);
        FFTW.fftwExecute(plan4);
        for (m in roundtrip3d) {
            Console.OUT.println("roundtrip" + m + " =" + roundtrip3d(m) / (K*K*K));
        }
    }
}
