/*
 * This file is part of ANU Molecular Mechanics (ANUMM).
 * ANUMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUMM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class DiatomicForceField implements ForceField {
    val diatomicPotentials : ValRail[DiatomicPotential];

    public def this(diatomicPotentials : ValRail[DiatomicPotential]) {
        this.diatomicPotentials = diatomicPotentials;
    }
    
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double {
        var V : Double = 0.0;        
        for((p) in 0..diatomicPotentials.length-1) {
            val potential = diatomicPotentials(p);
            V += potential.getPotentialAndForces();
        }
        return V;
    }
}
