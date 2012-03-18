/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2012.
 */
package au.edu.anu.mm;

import x10.util.Random;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

/**
 * Simulates an ion packet travelling in a circular path in a Penning trap
 * for the purposes of Fourier Transform Ion Cyclotron Resonance (FT-ICR) mass
 * spectroscopy.
 * Parameters to match simulation reported in:
 * @see Seung-Jin Han and Seung Koo Shin (1997)
 *  "Space-Charge Effects and Fourier transform ion cyclotron resonance signals:
 *   experimental observations and three-dimensional trajectory simulations"
 *   J. Am. Soc. Mass Spectrometry 8 (4), 319-326 
 * @author milthorpe
 */
public class TestCyclotron {
    public static def main(args : Array[String](1)) {
        var B:Double = 0.7646;
        var dt:Double = 50.0; // timestep in ns
        var V:Double = 1.0;
        var timesteps:Int = 80000; // number of timesteps
        var logSteps:Int = 100;
        var N:Int=3;
        if (args.size > 0) {
            B = Double.parseDouble(args(0));
            if (args.size > 1) {
                dt = Double.parseDouble(args(1));
                if (args.size > 2) {
                    V = Double.parseDouble(args(2));
                    if (args.size > 3) {
                        timesteps = Int.parseInt(args(3));
                        if (args.size > 4) {
                            logSteps = Int.parseInt(args(4));
                            if (args.size > 5) {
                                N = Int.parseInt(args(5));
                            }
                        }
                    }
                }
            }
        }

        Console.OUT.println("# Testing cyclotron: trapping potential: " + V + " magnetic field: " + B);


        val v = 3.4; // aiming for r ~ 2mm by r = mv/|q|B
        val q = 1.0;

        val species1 = "CH3CO";
        val mass1 = 43.04462;
        val f1 = q * B / mass1 * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r1 = mass1 * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %s predicted f = %10.2f r = %8.2f nm\n", species1, f1, r1*1e9);

        val species2 = "HCO";
        val mass2 = 29.0182;
        val f2 = q * B / mass2 * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r2 = mass2 * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %s predicted f = %10.2f r = %8.2f nm\n", species2, f2, r2*1e9);

        val rand = new Random();


        // initial distribution for each species is uniform cylinder along z dimension
        // centred at calculated radius for given velocity
        val atoms = new Array[MMAtom](N);
        for (i in 0..(N-1)) {
            val r = perturbation(rand);
            val theta = rand.nextDouble() * Math.PI * 2.0;
            val ex = Math.cos(theta) * r;
            val ey = Math.sin(theta) * r;
            val ion:MMAtom;
            if (i % 3 == 0) {
                ion = new MMAtom(species1, Point3d(-r1+ex, ey, perturbation(rand)), mass1, q);
            } else {
                ion = new MMAtom(species2, Point3d(-r2+ex, ey, perturbation(rand)), mass2, q);
            }

            // dominant velocity in y direction, slightly perturbed in random direction in X-Y plane
            // random velocity -1/2..1/2 in z direction
            val phi = rand.nextDouble() * Math.PI * 2.0;
            val evx = Math.cos(phi) * 1.0e-1 - 5.0e-2;
            val evy = Math.sin(phi) * 1.0e-1 - 5.0e-2;
            ion.velocity = Vector3d(evx, v+evy, rand.nextDouble()-0.5);
            atoms(i) = ion;
        }

        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val trap = new PenningTrap(N, distAtoms, V, new Vector3d(0.0, 0.0, B));
        trap.mdRun(dt, timesteps, logSteps);
    }

    private static def perturbation(rand:Random) {
        return rand.nextDouble()*1.0e-6 - 5.0e-7;
    }
}

