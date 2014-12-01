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
    public def getAtomTypes():Rail[AtomType] {
        throw new UnsupportedOperationException("GenericForceField.getAtomTypes");
    }

    public def getAtomMass(species:Int):Double {
        throw new UnsupportedOperationException("GenericForceField.getAtomMass");
    }

    public def getBondTypeIndex(species1:Int, species2:Int, structureFileType:Int):Int {
        throw new UnsupportedOperationException("GenericForceField.getBondTypeIndex");
    }
    
    public def getPotentialAndForces(particleDataPlh:PlaceLocalHandle[ParticleData]):Double {
        val energy = finish(Reducible.SumReducer[Double]()) {
            ateach(p in Dist.makeUnique()) {
                val particleData = particleDataPlh();

                var myEnergy:Double = 0.0;

                myEnergy += bondStretching(particleData);
                myEnergy += bondAngles(particleData);

                // TODO torsion, inversion, non-bonded terms
                offer myEnergy;
            }
        };
        return energy;
    }

    /** 
     * Compute all bond-stretch interactions between pairs of atoms joined by
     * covalent bonds.
     * TODO interaction potentials other than harmonic
     */
    private def bondStretching(particleData:ParticleData):Double {
        var potential:Double = 0.0;
        val bonds = particleData.bonds;
        if (bonds != null) {
            for (bond in bonds) {
                val bondParams = particleData.bondTypes(bond.typeIndex);
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
        }

        return potential;
    }

    /**
     * Compute all bond angle deformation interactions between triplets of
     * atoms, two of which share bonds with a common central atom.
     * The angle-bend interaction is described by a harmonic potential
     * proportional to the squared deviation of the angle from the 'ideal' bond
     * angle in radians: V = k (angle - ideal)^2.
     * The restorative force on the central atom is towards the interior of the
     * angle, and twice the strength of the force on each of the outside atoms.
     * TODO check direction of force for outside atoms
     * TODO interaction potentials other than harmonic
     */
    private def bondAngles(particleData:ParticleData):Double {
        var potential:Double = 0.0;
        val angles = particleData.angles;
        if (angles != null) {
            for (angle in angles) {
                val angleParams = particleData.angleTypes(angle.typeIndex);
                val atom1Center = particleData.x(angle.atom1Index);
                val atom2Center = particleData.x(angle.atom2Index);
                val atom3Center = particleData.x(angle.atom3Index);
                //Console.OUT.println("atom1 " + atom1Center + " atom2 " + atom2Center + " atom3 " + atom3Center);

                val r12 = atom2Center - atom1Center;
                val r32 = atom2Center - atom3Center;
                val dotProd = r12.dot(r32);
                val norm12 = r12.length();
                val norm32 = r32.length();
                val theta = Math.acos(dotProd / (norm12 * norm32));
                //Console.OUT.printf("ideal angle %10.6f theta %10.6f force constant %10.6f\n", angleParams.idealAngle, theta, angleParams.forceConstant);
                val bisectVector = (r12 + r32).normalize();
                val bend = theta - angleParams.idealAngle;
                val force = bend * angleParams.forceConstant * bisectVector;
                //Console.OUT.println("force = " + force);
                particleData.fx(angle.atom1Index) -= force;
                particleData.fx(angle.atom3Index) -= force;
                particleData.fx(angle.atom2Index) += 2.0 * force;

                // potential is for both particles
                val v = angleParams.forceConstant * bend * bend;
                potential += v;
            }
        }

        return potential;
    }
}
