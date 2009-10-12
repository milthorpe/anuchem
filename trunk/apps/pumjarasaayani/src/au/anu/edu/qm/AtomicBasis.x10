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
    global var orbitals:ArrayList[Orbital{self.at(this)}]{self.at(this)};

    public def this() { }

    public def make() { 
       orbitals = new ArrayList[Orbital{self.at(this)}]();
    } 

    public def addOrbital(orb:Orbital{self.at(this)}) {
       orbitals.add(orb);
    }

    public def getOrbitals() : ArrayList[Orbital{self.at(this)}]{self.at(this)} = this.orbitals;
}

