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
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

import x10.util.*;

/**
 * AtomicBasis.x10
 *
 * Represents an atomic basis as a collection of Orbital objects
 *
 * @author: V.Ganesh
 */
public class AtomicBasis { 
    global val orbitals:ArrayList[Orbital{self.at(this)}]{self.at(this)};

    public def this() { 
       orbitals = new ArrayList[Orbital{self.at(this)}]();
    } 

    public def addOrbital(orb:Orbital{self.at(this)}) {
       orbitals.add(orb);
    }

    public def getOrbitals() : ArrayList[Orbital{self.at(this)}]{self.at(this)} = this.orbitals;
}

