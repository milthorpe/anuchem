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

import x10.util.HashMap;
import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.mm.ForceField;
import au.edu.anu.mm.LennardJonesParameters;
import au.edu.anu.chem.BondType;
import au.edu.anu.chem.mm.MMAtom;

public class UniversalForceField implements ForceField {
    val BOND_ORDER_PROPORTIONALITY_CONSTANT = -0.1332;

    val defaultParams = UffParameters("default", 1.0, 0.0, 0.0, LennardJonesParameters("default", 0.0, 0.0, 0.0), 0.0, 0.0);
    val atomParameters : HashMap[String, UffParameters];

    /* The potential of the system as calculated by this force field.
       TODO should be shared var within getPotentialAndForces
     */
    var energy : Double = 0.0;

    public def this() {
        atomParameters = new HashMap[String, UffParameters]();
        // TODO read in parameters from file
        val hydrogen = UffParameters("H", 1.00794, 0.354, 180.0, LennardJonesParameters("H_", 12.0, 2.886, 0.044), 0.712, 4.53);
        atomParameters.put("H", hydrogen);
        val oxygen3 = UffParameters("O", 15.9994, 0.658, 104.51, LennardJonesParameters("O3", 14.085, 3.500, 0.060), 2.300, 8.741);
        atomParameters.put("O", oxygen3);
    }

    public def getAtomMass(symbol : String) : Double {
        return atomParameters.getOrElse(symbol, defaultParams).mass;
    }
    
    public def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double {
        val energy = finish(SumReducer()) {
            ateach(p in atoms) { 
                var myEnergy : Double = 0.0;
                val myAtoms = atoms(p);
                for([i] in 0..myAtoms.length()-1) {
                    val atomI = myAtoms(i);
                    atomI.force = Vector3d.NULL;
                    // bond stretching
                    //Console.OUT.println("atom " + atomI);
                    for (bond in atomI.getBonds()) {
                        if (bond.first.isStrongBond()) {
                            //Console.OUT.println("found bond: " + bond);
                            val atomJ = bond.second as MMAtom;
                            val paramsI = atomParameters.getOrElse(atomI.symbol, defaultParams);
                            val paramsJ = atomParameters.getOrElse(atomJ.symbol, defaultParams);
                            val bondStretch = getBondStretchTerm(bond.first, atomI, paramsI, atomJ, paramsJ);
                            myEnergy += bondStretch;
                        }
                    }
                }
                offer myEnergy;
            }
        };
        return energy;
    }

    /**
     * @param bond the bond type e.g. single, double, aromatic
     * @param atomI first atom
     * @param atomJ second atom
     * @param paramsI the UFF parameters for the first atom
     * @param paramsI the UFF parameters for the second atom
     * @return the bond stretch contribution (in Hartrees)
     */
    private def getBondStretchTerm(bond : BondType,
                                               atomI : MMAtom, paramsI : UffParameters,
                                               atomJ : MMAtom, paramsJ : UffParameters) {
        val direction = atomJ.centre - atomI.centre;
        val distance = direction.lengthSquared();
        val naturalDistance = getNaturalBondRadius(bond, atomI, paramsI, atomJ, paramsJ);
        //Console.OUT.println("naturalDistance = " + naturalDistance);

        // eq. 6 of UFF
        val forceConstant = paramsI.effectiveCharge * paramsJ.effectiveCharge /
                (naturalDistance * naturalDistance * naturalDistance);
        //Console.OUT.println("direction: " + direction + " forceConstant = " + forceConstant);
        val stretch = distance - naturalDistance;

        val force = direction.normalize() * (2.0 * forceConstant * stretch);
        //Console.OUT.println("bond stretching force: " + force);
        atomI.force = atomI.force + force;
        atomJ.force = atomJ.force + force.negate();

        // eq. 1a of UFF
        return forceConstant * stretch * stretch;
    }

    /**
     * Get the natural bond radius for bonded atoms, which is the sum of the
     * bond radii of the two atoms, plus a bond order correction, minus an
     * electronegativity correction.  See eq. 2 of UFF.
     * NOTE: there is a typo in the UFF paper, it suggests the electronegativity
     * correction should be added.
     * see http://towhee.sourceforge.net/forcefields/uff.html
     * @param bond the bond type e.g. single, double, aromatic
     * @param atomI first atom
     * @param atomJ second atom
     * @param paramsI the UFF parameters for the first atom
     * @param paramsI the UFF parameters for the second atom
     * @return the natural bond radius between the atoms
     */
    private def getNaturalBondRadius(bond : BondType, 
                                                atomI : MMAtom, paramsI : UffParameters,
                                               atomJ : MMAtom, paramsJ : UffParameters) {
        val radiusI = paramsI.bondRadius;
        val radiusJ = paramsJ.bondRadius;
        val naturalRadius = radiusI + radiusJ;
        //Console.OUT.println("naturalRadius = " + naturalRadius);
        // eq. 3 of UFF
        val bondOrderCorrection = BOND_ORDER_PROPORTIONALITY_CONSTANT * naturalRadius * Math.log(bond.bondOrder);
        //Console.OUT.println("bondOrderCorrection = " + bondOrderCorrection);
        val chiI = paramsI.electronegativity;
        val chiJ = paramsJ.electronegativity;
        val dSqrtChi = Math.sqrt(chiI) - Math.sqrt(chiJ);
        // eq. 4 of UFF
        val electronegativityCorrection = radiusI * radiusJ * (dSqrtChi * dSqrtChi) / (chiI * radiusI + chiJ * radiusJ);

        return (naturalRadius + bondOrderCorrection - electronegativityCorrection);
    }

    static struct SumReducer implements Reducible[Double] {
        public def zero() = 0.0;
        public def apply(a:Double, b:Double) = (a + b);
    }
}
