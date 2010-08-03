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
 * Tests distributed Fast Fourier Transform.
 */
public class TestDistributedFFT {
    public static def main(args : Rail[String]!) {
        testFft3d();
    }

    public static def testFft3d() {
        val N = 3;

        val twoPlusI = 2.0 + Complex.I;
        val r = Region.make(0, N-1);
        val r3 = (r * r * r) as Region(3);
        val input = DistArray.make[Complex](Dist.makeBlockBlock(r3, 0, 1), (p(i,j,k) : Point(3)) => twoPlusI * (1.0 / (i + j + k + 1)));
        Console.OUT.println("dist = " + input.dist);
        val output = DistArray.make[Complex](input.dist);
        val temp = DistArray.make[Complex](input.dist);
        Console.OUT.println("input");
        for (p in input.dist.places()) {
            at (p) {
                for (p in input | here) {
                    Console.OUT.println(input(p));
                }
            }
        }

        val fft = new Distributed3dFft(N);

        //Console.OUT.println("\noutput");
        fft.doFFT3d(input, output, temp, false);
        for (p in output.dist.places()) {
            at (p) {
                for (p in output | here) {
                    //Console.OUT.println(output(p));
                }
            }
        }

        Console.OUT.println("\nroundtrip");
        fft.doFFT3d(output, output, temp, true);
        for (p in output.dist.places()) {
            at (p) {
                for (p in output | here) {
                    Console.OUT.println((output(p) / (N*N*N)));
                }
            }
        }
    }
}
