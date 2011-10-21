/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/lipcenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.fft;

/**
 * Tests distributed Fast Fourier Transform.
 */
public class TestDistributed3dFft {
    public static def main(args : Array[String](1)) {
        testFft3d();
    }

    public static def testFft3d() {
        val N = 32;

        val twoPlusI = 2.0 + Complex.I;
        val r = Region.make(0, N-1);
        val r3 = (r * r * r) as Region(3);
        val input = DistArray.make[Complex](Dist.makeBlockBlock(r3, 0, 1), (p[i,j,k] : Point(3)) => twoPlusI * (1.0 / (i + j + k + 1)));
        Console.OUT.println("dist = " + input.dist);
        val output = DistArray.make[Complex](input.dist);
        val temp = DistArray.make[Complex](input.dist);
        Console.OUT.println("input");
        for (p1 in input.dist.places()) at(p1) {
            for (p in input.dist | here) {
                //Console.OUT.println(input(p));
            }
        }

        val fft = new Distributed3dFft(N, input, output, temp);

        //Console.OUT.println("\noutput");
        fft.doFFT3d(false);
        for (p1 in output.dist.places()) at(p1) {
            for (p in output.dist | here) {
                //Console.OUT.println(output(p));
            }
        }

        val roundtrip = new Distributed3dFft(N, output, output, temp);

        Console.OUT.println("\nroundtrip");
        roundtrip.doFFT3d(true);
        for (p1 in output.dist.places()) at(p1) {
            for (p in output.dist | here) {
                //Console.OUT.println((output(p) / (N*N*N)));
                if (Math.abs((output(p)/(N*N*N)).re - input(p).re) > 1.0e-10) {
                    Console.OUT.println("input was " + input(p));
                    Console.OUT.println((output(p)/(N*N*N)).re - input(p).re);
                }
            }
        }
    }
}