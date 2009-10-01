/**
 * AtomicBasis.x10
 *
 * Represents an atomic basis as a collection of Orbital objects
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class AtomicBasis { 
    val orbitals:ArrayList[Orbital];

    public def this() { 
       orbitals = new ArrayList[Orbital]();
    } 

    public def addOrbital(orb:Orbital) {
       orbitals.add(orb);
    }

    public def getOrbitals() = this.orbitals;
}

