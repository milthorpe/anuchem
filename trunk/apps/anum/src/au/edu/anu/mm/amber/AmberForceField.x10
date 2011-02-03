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
package au.edu.anu.mm.amber;

import x10.util.HashMap;
import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.mm.ForceField;
import au.edu.anu.chem.mm.MMAtom;

public class AmberForceField implements ForceField {
    val bondParameters : HashMap[String, AmberBondParameters];
    val angleParameters : HashMap[String, AmberAngleParameters];

    /* The potential of the system as calculated by this force field.
       TODO should be shared var within getPotentialAndForces
     */
    var energy : Double = 0.0;

    public def this() {
        bondParameters = new HashMap[String, AmberBondParameters]();
        // TODO read in parameters from file
        val waterBond = AmberBondParameters("OW-HW", 553.0, 0.9572);
        bondParameters.put(waterBond.description, waterBond);

        angleParameters = new HashMap[String, AmberAngleParameters]();
        // TODO read in parameters from file
        val waterAngle = AmberAngleParameters("OW-HW", 100.0, 104.52);
        angleParameters.put(waterAngle.description, waterAngle);
    }
    
    public def getPotentialAndForces(atoms: DistArray[Array[MMAtom](1){rect,rail}](1)) : Double {
        energy = 0.0;
        finish ateach(p in atoms) { 
            var myEnergy : Double = 0.0;
            val myAtoms = atoms(p);
            for([i] in 0..myAtoms.length()-1) {
                val atom = myAtoms(i);
                atom.force = Vector3d.NULL;
                // bond stretching
                Console.OUT.println("atom " + atom);
                for (bond in atom.getBonds()) {
                    if (bond.first.isStrongBond()) {
                        Console.OUT.println("found bond: " + bond);
                        var bondDescription : String;
                        if (atom.symbol < bond.second.symbol) {
                            bondDescription = atom.symbol + "-" + bond.second.symbol;
                        } else {
                            bondDescription = bond.second.symbol + "-" + atom.symbol;
                        }
                        Console.OUT.println("type: " + bondDescription);
                        val bondParam = bondParameters.getOrElse(bondDescription, AmberBondParameters("null", 0.0, 0.0));
                        val r = atom.centre - bond.second.centre;
                        val bondRadius = Math.abs(r.length());
                        Console.OUT.println("radius = " + bondRadius + ", ideal = " + bondParam.bondRadius);
                        val displacement = bondRadius - bondParam.bondRadius;
                        val force = 2 * bondParam.forceConstant * displacement * r.normalize();
                        atom.force = atom.force + force;
                        Console.OUT.println("force = " + force);
                        val potential = bondParam.forceConstant * displacement * displacement;
                        Console.OUT.println("potential = " + potential);
                        myEnergy += potential;
                    }
                }
            }
            val myEnergyFinal = myEnergy;
            at (this) { atomic { energy += myEnergyFinal; } };
        }
        return energy;
    }
}
