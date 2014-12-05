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

import x10.util.ArrayList;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

import au.edu.anu.chem.mm.ParticleData;

public class AnummParticleData extends ParticleData {
    public var boxEdges:Vector3d;

    // Molecule / residue data

    /** The different molecule types found in the simulation system */
    public var moleculeTypes:Rail[MoleculeType];
    /** Mapping from residue number to molecule type */
    public var moleculeTypeIndex:Rail[Int];

    /** Residue number (e.g. index of molecule) for each atom */
    public val residueNumber = new ArrayList[Int]();
    /** Atom velocities */
    public val dx = new ArrayList[Vector3d]();

    // Bond data
    public var bondTypes:Rail[BondStretchParameters];
    public var angleTypes:Rail[BondAngleParameters];

    /** Bonds between atoms, where the first atom is held locally at this place */
    public var bonds:Rail[Bond];
    /** Bond angles between atoms, where the second atom is held locally at this place */
    public var angles:Rail[BondAngle];

    public def addAtom(index:Long, atomType:Int, center:Point3d) {
        super.addAtom(index, atomType, center);
        residueNumber.add(0n);
        dx.add(Vector3d.NULL);
    }
}

