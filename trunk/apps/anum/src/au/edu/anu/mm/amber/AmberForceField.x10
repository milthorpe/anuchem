/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
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
    global val bondParameters : HashMap[String, AmberBondParameters];
    global val angleParameters : HashMap[String, AmberAngleParameters];

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
    
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double {
        energy = 0.0;
        finish ateach(p in atoms) { 
            var myEnergy : Double = 0.0;
            val myAtoms = atoms(p);
            for((i) in 0..myAtoms.length()-1) {
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
