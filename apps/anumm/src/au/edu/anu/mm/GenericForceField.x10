/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10.regionarray.Dist;

import au.edu.anu.chem.mm.MMAtom;

public class GenericForceField implements ForceField {
    val specs:Rail[AtomType];
    val bondStretchParameters:Rail[BondStretchParameters];

    public def this() {
        // TODO read in parameters from file
        specs = new Rail[AtomType](3);
        specs(0) = AtomType("OW", 8n, 15.9949146, -0.82);
        specs(1) = AtomType("HW1", 1n, 1.00794,    0.41);
        specs(2) = AtomType("HW2", 1n, 1.00794,    0.41);

        bondStretchParameters = new Rail[BondStretchParameters](2);
        bondStretchParameters(1) = BondStretchParameters(0n, 1n, 0.1, 418400.0);
    }

    public def getAtomTypes() {
        return specs;
    }

    public def getAtomMass(species:Int):Double {
        return specs(species).mass;
    }

    public def getBondTypeIndex(species1:Int, species2:Int, structureFileType:Int):Int {
        // TODO lookup bond stretch params
        return 0n;
    }
    
    public def getPotentialAndForces(particleDataPlh:PlaceLocalHandle[ParticleData]):Double {
        val energy = finish(SumReducer()) {
            ateach(p in Dist.makeUnique()) {
                val particleData = particleDataPlh();

                var myEnergy : Double = 0.0;

                myEnergy += bondStretching(particleData);

                // TODO angle, torsion, inversion, non-bonded terms
                offer myEnergy;
            }
        };
        return energy;
    }

    private def bondStretching(particleData:ParticleData):Double {
        var potential:Double = 0.0;
        val bonds = particleData.bonds;
        for (bond in bonds) {
            val bondParams = bondStretchParameters(bond.typeIndex);
            val atom1Center = particleData.x(bond.atom1Index);
            val atom2Center = particleData.x(bond.atom2Index);

            //Console.OUT.println("atom1 species " + bondParams.species1 + " " + atom1Center);
            //Console.OUT.println("atom2 species " + bondParams.species2 + " " + atom2Center);
            //Console.OUT.printf("ideal radius %10.6f forceConstant %10.6f\n", bondParams.idealRadius, bondParams.forceConstant);

            val r12 = atom2Center - atom1Center;
            val bondLength = r12.length();

            val stretch = bondLength - bondParams.idealRadius;
            //Console.OUT.printf("stretch %10.6f r12 %10.6f %10.6f %10.6f\n", stretch, r12.i, r12.j, r12.k);
            val normalizedVec = r12 * (1.0 / bondLength);
            val force = normalizedVec * (bondParams.forceConstant * stretch);
            //Console.OUT.printf("bond stretching force %10.6f %10.6f %10.6f\n", force.i, force.j, force.k);

            particleData.fx(bond.atom1Index) += force;
            particleData.fx(bond.atom2Index) -= force;
            //Console.OUT.println("atom1 force " + particleData.fx(bond.atom1Index));
            //Console.OUT.println("atom2 force " + particleData.fx(bond.atom2Index));

            // potential is for both particles
            val v = bondParams.forceConstant * stretch * stretch;
            //Console.OUT.println("potential = " + v);
            potential += v;
        }

        return potential;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
