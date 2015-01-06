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

import x10.util.Pair;

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

    // Bond data
    public var bondTypes:Rail[BondStretchParameters];
    public var angleTypes:Rail[BondAngleParameters];

    /** Bonds between atoms, where the first atom is held locally at this place */
    public var bonds:Rail[Bond];
    /** Bond angles between atoms, where the second atom is held locally at this place */
    public var angles:Rail[BondAngle];

    /** Get the maximum distance from the origin in any dimension for any atom. */
    public def getMaxExtent():Vector3d {
        var maxX:Double = 0.0;
        var maxY:Double = 0.0;
        var maxZ:Double = 0.0;
        for (center in x) {
            maxX = Math.max(maxX, Math.abs(center.i));
            maxY = Math.max(maxY, Math.abs(center.j));
            maxZ = Math.max(maxX, Math.abs(center.k));
        }
        return Vector3d(maxX, maxY, maxZ);
    }

    /**
     * Create the list of exclusions for non-bonded interactions.
     * All atoms that are separated by at most three bonds are excluded
     * from non-bonded interactions.
     */
    public def generateNonBondedExclusions() {
        if (bonds != null) {
            for (bond in bonds) {
                exclusions.add(Pair[Long,Long](
                    globalIndex(bond.atom1Index),
                    globalIndex(bond.atom2Index)));
                exclusions.add(Pair[Long,Long](
                    globalIndex(bond.atom2Index),
                    globalIndex(bond.atom1Index)));
            }
        }

        if (angles != null) {
            for (angle in angles) {
                exclusions.add(Pair[Long,Long](
                    globalIndex(angle.atom1Index),
                    globalIndex(angle.atom3Index)));
                exclusions.add(Pair[Long,Long](
                    globalIndex(angle.atom3Index),
                    globalIndex(angle.atom1Index)));
            }
        }

        exclusions.sort(compareExclusions);
    }
}

