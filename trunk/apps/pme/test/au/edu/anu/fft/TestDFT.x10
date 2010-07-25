/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
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
