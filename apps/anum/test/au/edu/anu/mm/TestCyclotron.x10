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
import au.edu.anu.mm.SystemProperties;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.TestElectrostatic;

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
public class TestCyclotron extends TestElectrostatic {
    public static def main(args : Array[String](1)) {
        var B:Double = 0.7646;
        var dt:Double = 50000.0; // timestep in fs
        var V:Double = 1.0;
        var timesteps:Int = 80000; // number of timesteps
        var logSteps:Int = 100;
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
                        }
                    }
                }
            }
        }

        Console.OUT.println("# Testing cyclotron: trapping potential: " + V + " magnetic field: " + B);


        val v = 3.4; // aiming for r ~ 2mm by r = mv/|q|B
        val q = 1.0;

        val species1 = "CH3CO";
        val f1 = q * B / PenningTrap.getAtomMass(species1) * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r1 = PenningTrap.getAtomMass(species1) * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %s predicted f = %10.2f r = %8.2f nm\n", species1, f1, r1*1e9);

        val ion1 = new MMAtom(species1, Point3d(-r1, 0.0, 0.0), q);
        ion1.velocity = Vector3d(0.0, v, -0.1);

        val species2 = "HCO";
        val f2 = q * B / PenningTrap.getAtomMass(species2) * (PenningTrap.CHARGE_MASS_FACTOR) / (2.0 * Math.PI);
        val r2 = PenningTrap.getAtomMass(species2) * v / (q * B) / PenningTrap.CHARGE_MASS_FACTOR;
        Console.OUT.printf("# %s predicted f = %10.2f r = %8.2f nm\n", species2, f2, r2*1e9);

        val ion2 = new MMAtom(species2, Point3d(-r2, 0.0, 0.0), q);
        ion2.velocity = Vector3d(0.0, v, -0.1);


        val atoms = new Array[MMAtom](2);
        atoms(0) = ion1;
        atoms(1) = ion2;

        val distAtoms = DistArray.make[Rail[MMAtom]](Dist.makeBlock(0..0, 0));
        distAtoms(0) = atoms;

        val trap = new PenningTrap(1, distAtoms, V, new Vector3d(0.0, 0.0, B));
        val kineticEnergy = (a:MMAtom) => PenningTrap.getAtomMass(a.symbol) * PenningTrap.MASS_FACTOR * a.velocity.lengthSquared() * 1e24;
        val potentialEnergy = (a:MMAtom) => a.charge*PenningTrap.CHARGE_FACTOR*trap.getElectrostaticPotential(a.centre) * 1e24;
        val totalEnergy = (a:MMAtom) => PenningTrap.getAtomMass(a.symbol) * PenningTrap.MASS_FACTOR * a.velocity.lengthSquared() * 1e24
                                        + a.charge * PenningTrap.CHARGE_FACTOR * trap.getElectrostaticPotential(a.centre) * 1e24;
        val detectorCurrent = (a:MMAtom) => PenningTrap.getImageCurrent(a) * PenningTrap.CHARGE_FACTOR * 1e18;

        val onePFs = [OneParticleFunction("mean_X (mm)", (a:MMAtom) => a.centre.i*1e6), 
                      OneParticleFunction("mean_Y (mm)", (a:MMAtom) => a.centre.j*1e6), 
                      OneParticleFunction("mean_Z (mm)", (a:MMAtom) => a.centre.k*1e6), 
                      OneParticleFunction("Ek (yJ)", kineticEnergy), 
                      OneParticleFunction("Ep (yJ)", potentialEnergy), 
                      OneParticleFunction("E (yJ)", totalEnergy),
                      OneParticleFunction("I (aA)", detectorCurrent)];
        trap.setProperties(new SystemProperties(atoms.size, onePFs));
        trap.mdRun(dt, timesteps, logSteps);
    }
}

