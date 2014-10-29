/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

public class ParticleData {
    public var description:String;

    // Atom data
    public var numAtoms:Long;
    /** Species of each atom */
    public var species:Rail[Int];
    /** Global index of each atom */
    public var globalIndex:Rail[Long];
    /** Atom centers in nm */
    public var x:Rail[Point3d];
    /** Atom velocities */
    public var dx:Rail[Vector3d];
    /** Instantaneous force on atoms */
    public var fx:Rail[Vector3d];

    // Bond data

    /** Bonds between atoms, where the first atom is held locally at this place */
    public var bonds:Rail[Bond];

    public def allocateAtoms(N:Long) {
        this.numAtoms = N;
        species = new Rail[Int](N);
        globalIndex = new Rail[Long](N);
        x = new Rail[Point3d](N);
        dx = new Rail[Vector3d](N);
        fx = new Rail[Vector3d](N);
    }
}

