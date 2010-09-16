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
 * Tests regular (non-fast) Discrete Fourier Transforms.
 */
public class TestDFT {
    public static def main(args : Array[String](1)) {
        testDft1D();
        testDft3D();
    }

    public static def testDft1D() {
        val N = 10;

        val onePlusI = Complex.ONE + Complex.I;
        val input = new Array[Complex](N, ([i] : Point) => onePlusI * (1.0 / (i + 1)));
        Console.OUT.println("input");
        for (var i:int=0; i<input.size; i++) {
            Console.OUT.println(input(i));
        }

        Console.OUT.println("\noutput");
        val output = DFT.dft1D(input, false);
        for (var i:int=0; i<output.size; i++) {
            Console.OUT.println(output(i));
        }

        Console.OUT.println("\nroundtrip");
        val roundtrip = DFT.dft1D(output, true);
        for (var i:int=0; i<roundtrip.size; i++) {
            Console.OUT.println(roundtrip(i) / N);
        }
    }

    public static def testDft3D() {
        val N = 3;

        val twoPlusI = 2.0 + Complex.I;
        val r = Region.make(0, N-1);
        val r3 = (r * r * r) as Region(3){rect,zeroBased};
        val input = new Array[Complex](r3, (p[i,j,k] : Point(3)) => twoPlusI * (1.0 / (i + j + k + 1)));
        Console.OUT.println("input");
        for (p in input) {
            Console.OUT.println(input(p));
        }

        Console.OUT.println("\noutput");
        val output = DFT.dft3D(input, false);
        for (p in output) {
            Console.OUT.println(output(p));
        }

        Console.OUT.println("\nroundtrip");
        val roundtrip = DFT.dft3D(output, true);
        for (p in roundtrip) {
            Console.OUT.println(roundtrip(p) / (N*N*N));
        }
    }
}
