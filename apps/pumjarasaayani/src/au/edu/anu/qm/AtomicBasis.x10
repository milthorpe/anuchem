/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2011.
 */
package au.edu.anu.qm;

/**
 * Represents an atomic basis as a collection of Orbital objects
 *
 * @author: V.Ganesh
 */
public class AtomicBasis { 
    public val orbitals:Rail[Orbital];

    public def this(orbitals:Rail[Orbital]) { 
       this.orbitals = orbitals;
    } 
}

