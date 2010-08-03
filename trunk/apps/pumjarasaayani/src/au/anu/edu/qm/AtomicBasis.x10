/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
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

