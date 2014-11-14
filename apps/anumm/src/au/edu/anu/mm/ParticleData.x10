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

    public var boxEdges:Vector3d;

    // Molecule / residue data

    /** The different molecule types found in the simulation system */
    public var moleculeTypes:Rail[MoleculeType];
    /** Mapping from residue number to molecule type */
    public var moleculeTypeIndex:Rail[Int];

    // Atom data

    public var atomTypes:Rail[AtomType];

    public var numAtoms:Long;
    /** Atom type index for each atom */
    public var atomTypeIndex:Rail[Int];
    /** Residue number (e.g. index of molecule) for each atom */
    public var residueNumber:Rail[Int];
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
        atomTypeIndex = new Rail[Int](N);
        residueNumber = new Rail[Int](N);
        globalIndex = new Rail[Long](N);
        x = new Rail[Point3d](N);
        dx = new Rail[Vector3d](N);
        fx = new Rail[Vector3d](N);
    }
}

