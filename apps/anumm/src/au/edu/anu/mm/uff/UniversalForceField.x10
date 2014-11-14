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
package au.edu.anu.mm.uff;

import x10.regionarray.Dist;
import x10.util.HashMap;
import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.mm.BondStretchParameters;
import au.edu.anu.mm.ForceField;
import au.edu.anu.mm.LennardJonesParameters;
import au.edu.anu.mm.ParticleData;
import au.edu.anu.mm.AtomType;
import au.edu.anu.chem.BondType;

public class UniversalForceField implements ForceField {
    public static SPECIES_H = 1;
    public static SPECIES_C = 6;
    public static SPECIES_N = 7;
    public static SPECIES_O = 8;

    val BOND_ORDER_PROPORTIONALITY_CONSTANT = -0.5573; // kJ/mol

    val atomParameters:Rail[UffParameters];

    val bondTypes:Rail[BondType] = [BondType.NO_BOND, 
                                    BondType.WEAK_BOND,
                                    BondType.SINGLE_BOND,
                                    BondType.DOUBLE_BOND,
                                    BondType.TRIPLE_BOND,
                                    BondType.QUADRUPLE_BOND,
                                    BondType.AROMATIC_BOND,
                                    BondType.AMIDE_BOND,
                                    BondType.IONIC_BOND];

    public static BOND_TYPE_INDEX_SINGLE   = 2n;
    public static BOND_TYPE_INDEX_DOUBLE   = 3n;
    public static BOND_TYPE_INDEX_TRIPLE   = 4n;
    public static BOND_TYPE_INDEX_AROMATIC = 6n;

    public def this() {
        atomParameters = new Rail[UffParameters](9);
        // TODO read in parameters from file
        atomParameters(SPECIES_H) = UffParameters("H", 1.00794, 3.54, 180.0,  LennardJonesParameters("H_",  12.0,   28.86, 0.044), 0.712, 4.53);
        atomParameters(SPECIES_C) = UffParameters("C", 12.0000, 7.57, 109.47, LennardJonesParameters("C_3", 12.73,  38.51, 0.105), 1.912, 5.34);
        atomParameters(SPECIES_N) = UffParameters("N", 14.0031, 7.00, 106.7,  LennardJonesParameters("N_3", 13.407, 36.60, 0.069), 2.544, 6.899);
        atomParameters(SPECIES_O) = UffParameters("O", 15.9994, 6.58, 104.51, LennardJonesParameters("O_3", 14.085, 35.00, 0.060), 2.300, 8.741);
    }

    public def getAtomTypes() {
        val specs = new Rail[AtomType](atomParameters.size);
        for (i in 0..(atomParameters.size-1)) {
            val atom = atomParameters(i);
            if (atom.description != null) {
                specs(i) = new AtomType(atom.description, i as Int, atom.mass, atom.effectiveCharge);
            }
        }
        return specs;
    }

    public def getAtomMass(species:Int) : Double {
        return atomParameters(species).mass;
    }

    public def getBondTypeIndex(species1:Int, species2:Int, structureFileType:Int):Int {
        var bondTypeIndex:Int;
        switch (structureFileType) {
            case (1n):
                bondTypeIndex = UniversalForceField.BOND_TYPE_INDEX_SINGLE;
                break;
            case (2n):
                bondTypeIndex = UniversalForceField.BOND_TYPE_INDEX_DOUBLE;
                break;
            case (3n):
                bondTypeIndex = UniversalForceField.BOND_TYPE_INDEX_TRIPLE;
                break;
            case (4n):
                bondTypeIndex = UniversalForceField.BOND_TYPE_INDEX_AROMATIC;
                break;
            default:
                Console.ERR.println("Bond type " + structureFileType + " not found. Defaulting to single bond.");
                bondTypeIndex = UniversalForceField.BOND_TYPE_INDEX_SINGLE;
        }
        return bondTypeIndex;
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
            val atom1Center = particleData.x(bond.atom1Index);
            val atom2Center = particleData.x(bond.atom2Index);

            val bondParams = getBondStretchParameters(particleData.atomTypeIndex(bond.atom1Index), particleData.atomTypeIndex(bond.atom2Index), bond.typeIndex);

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

    private def getBondStretchParameters(species1:Int, species2:Int, bondTypeIndex:Int) {
        val bondType = bondTypes(bondTypeIndex);
        val paramsI = atomParameters(species1);
        val paramsJ = atomParameters(species2);

        /*
         * Get the natural bond radius for bonded atoms, which is the sum of the
         * bond radii of the two atoms, plus a bond order correction, minus an
         * electronegativity correction.  See eq. 2 of UFF.
         * NOTE: there is a typo in the UFF paper, it suggests the electronegativity
         * correction should be added.
         * see http://towhee.sourceforge.net/forcefields/uff.html
         */
        val radiusI = paramsI.bondRadius;
        val radiusJ = paramsJ.bondRadius;
        val naturalRadius = radiusI + radiusJ;
        //Console.OUT.println("naturalRadius = " + naturalRadius);
        // eq. 3 of UFF
        val bondOrderCorrection = BOND_ORDER_PROPORTIONALITY_CONSTANT * naturalRadius * Math.log(bondType.bondOrder);
        //Console.OUT.println("bondOrderCorrection = " + bondOrderCorrection);
        val chiI = paramsI.electronegativity;
        val chiJ = paramsJ.electronegativity;
        val dSqrtChi = Math.sqrt(chiI) - Math.sqrt(chiJ);
        // eq. 4 of UFF
        val electronegativityCorrection = radiusI * radiusJ * (dSqrtChi * dSqrtChi) / (chiI * radiusI + chiJ * radiusJ);

        val idealRadius = (naturalRadius + bondOrderCorrection - electronegativityCorrection);

        // eq. 6 of UFF
        val forceConstant = paramsI.effectiveCharge * paramsJ.effectiveCharge /
                (idealRadius * idealRadius * idealRadius);

        return BondStretchParameters(species1, species2, idealRadius, forceConstant);
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public operator this(a:Double, b:Double) = (a + b);
    }
}
