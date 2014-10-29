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

import au.edu.anu.mm.ForceField;
import au.edu.anu.mm.SpeciesSpec;
import au.edu.anu.chem.mm.MMAtom;

public class GenericForceField implements ForceField {
    val specs:Rail[SpeciesSpec];
    val bondStretchParameters:Rail[BondStretchParameters];

    public def this() {
        // TODO read in parameters from file
        specs = new Rail[SpeciesSpec](3);
        specs(0) = SpeciesSpec("OW", 15.9949146, -0.82, 8n);
        specs(1) = SpeciesSpec("HW1", 1.00794,    0.41, 1n);
        specs(2) = SpeciesSpec("HW2", 1.00794,    0.41, 1n);

        bondStretchParameters = new Rail[BondStretchParameters](1);
        bondStretchParameters(0) = BondStretchParameters(0n, 1n, 0.1, 418400.0);
    }

    public def getSpeciesSpecs() {
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
        var energy:Double = 0.0;
        val bonds = particleData.bonds;
        for (bond in bonds) {
            val bondParams = bondStretchParameters(bond.typeIndex);
            val atom1Center = particleData.x(bond.atom1Index);
            val atom2Center = particleData.x(bond.atom2Index);

            //Console.OUT.println("atom1 species " + bondParams.species1 + " " + atom1Center);
            //Console.OUT.println("atom2 species " + bondParams.species2 + " " + atom2Center);
            //Console.OUT.printf("ideal radius %10.6f forceConstant %10.6f\n", bondParams.idealRadius, bondParams.forceConstant);

            val direction = atom2Center - atom1Center;
            val distance = direction.lengthSquared();

            val stretch = distance - bondParams.idealRadius;
            //Console.OUT.printf("stretch  %10.6f direction %10.6f %10.6f %10.6f\n", stretch, direction.i, direction.j, direction.k);
            val force = direction.normalize() * (bondParams.forceConstant * stretch);
            //Console.OUT.printf("bond stretching force %10.6f %10.6f %10.6f\n", force.i, force.j, force.k);

            particleData.fx(bond.atom1Index) += force;
            particleData.fx(bond.atom2Index) -= force;

            energy += 0.5 * bondParams.forceConstant * stretch * stretch;
        }

        return energy;
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
