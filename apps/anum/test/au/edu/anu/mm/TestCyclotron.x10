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
 * @see Leach et al. (2010). "Comparison of Particle-In-Cell simulations with 
 *   experimentally observed frequency shifts between ions of the same mass-to-
 *   charge in Fourier Transform Ion Cyclotron Resonance Mass Spectrometry".
 *   J. Am. Soc. Mass Spectrometry 21 (2), 203-208 
 * @author milthorpe
 */
public class TestCyclotron {
    public static def main(args : Array[String](1)) {
        var N:Int=100; // number of particles
        var timesteps:Int = 16000; // number of timesteps
        var dt:Double = 25.0; // timestep in ns
        var logSteps:Int = 1;
        var B:Double = 7.0; // magnetic field
        var V:Double = 1.0; // trapping potential
        if (args.size > 0) {
            N = Int.parseInt(args(0));
            if (args.size > 1) {
                timesteps = Int.parseInt(args(1));
                if (args.size > 2) {
                    dt = Double.parseDouble(args(2));
                    if (args.size > 3) {
                        logSteps = Int.parseInt(args(3));
                        if (args.size > 4) {
                            B = Double.parseDouble(args(4));
                            if (args.size > 5) {
                                V = Double.parseDouble(args(5));
                            }
                        }
                    }
                }
            }
        }

        Console.OUT.printf("# Testing cyclotron: trapping potential: %2.1f V magnetic field: %3.1f T\n", V, B);


        val v = 2.7e4; // aiming for r ~ 6mm by r = mv/|q|B
        val q = 1.0;

        val species1 = "Cs+"; // "CH3CO";
        val mass1 = 132.9054; // 43.04462;
        val f1 = q * B / mass1 * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r1 = mass1 * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %8s predicted f = %8i Hz r = %6.3f mm\n", species1, f1 as Int, r1*1e3);

        val species2 = "x150"; //"HCO";
        val mass2 = 150.00; // 29.0182;
        val f2 = q * B / mass2 * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r2 = mass2 * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %8s predicted f = %8i Hz r = %6.3f mm\n", species2, f2 as Int, r2*1e3);

        val rand = new Random(27178281L);

        // initial distribution for each species is uniform 1mm cylinder along z dimension
        // centred at calculated radius for given velocity
        val atoms = new Array[MMAtom](N);
        for (i in 0..(N-1)) {
            val er = rand.nextDouble() * 1.0e-3;
            val theta = rand.nextDouble() * Math.PI * 2.0;
            val ex = Math.cos(theta) * er;
            val ey = Math.sin(theta) * er;
            val ion:MMAtom;
            if (i % 2 == 0) {
                ion = new MMAtom(species1, Point3d(-r1+ex, ey, perturbation(rand, 1e-3)), mass1, q);
            } else {
                ion = new MMAtom(species2, Point3d(-r2+ex, ey, perturbation(rand, 1e-3)), mass2, q);
            }

            // velocity is Maxwellian distribution with addition of velocity v in y direction
            val ev = maxwellianVelocity(rand, ion.mass);
            ion.velocity = Vector3d(ev.i, v+ev.j, ev.k);
            atoms(i) = ion;
        }

        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val trap = new PenningTrap(N, distAtoms, V, new Vector3d(0.0, 0.0, B));
        trap.mdRun(dt, timesteps, logSteps);
    }

    private static def perturbation(rand:Random, max:Double) {
        return rand.nextDouble()*max*2.0 - max;
    }

    /** 
     * Returns a velocity vector taken from a Maxwellian distribution at 300K,
     * i.e. a vector with each component a normal distribution with mean 0 and
     * variance kT/m
     * @param rand a uniform random number generator
     * @param mass the particle mass
     */
    private static def maxwellianVelocity(rand:Random, mass:Double):Vector3d {
        val variance = 8.314462175 * 300 / mass; // (gas constant R) * 300K / (molecular mass m)
        return Vector3d(randNormal(rand, variance), randNormal(rand, variance), randNormal(rand, variance));
    }

    /**
     * Use the Box-Muller transform to generate a pseudo-random number from a normal distribution.
     * @param rand a uniform random number generator
     * @param variance distribution variance 
     * @return a random variable taken from a normal distribution with mean 0
     */
    private static def randNormal(rand:Random, variance:Double) {
        val u1 = rand.nextDouble();
        val u2 = rand.nextDouble();
        val z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        return variance * z0;
    }
}

