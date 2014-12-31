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

public class MoleculeType {
    public var name:String;

    // TODO link atomName and index
    public var atomNames:ArrayList[String];
    public var atomTypes:ArrayList[Int];

    public var bonds:ArrayList[Bond];
    public var angles:ArrayList[BondAngle];
    // TODO angles, dihedrals
}

