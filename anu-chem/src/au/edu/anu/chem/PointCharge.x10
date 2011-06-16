/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011.
 */
package au.edu.anu.chem;

import x10x.vector.Point3d;

/**
 * This class represents a moveable point charge, which may
 * be used for the purposes of electrostatic calculation in
 * a molecular dynamics simulation.
 *
 * @author milthorpe
 */
public struct PointCharge { 
    public val centre : Point3d;
    public val charge : Double;

	public def this(centre : Point3d, charge : Double) {
		this.centre = centre;
		this.charge = charge;
	}

    public def toString() : String {
        return "Point charge " + charge + " " + centre;
    }
}

