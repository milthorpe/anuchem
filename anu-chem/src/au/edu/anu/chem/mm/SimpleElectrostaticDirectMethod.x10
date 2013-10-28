/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2013.
 */
package au.edu.anu.chem.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.vector.Tuple3d;
import au.edu.anu.util.Timer;
import au.edu.anu.chem.PointCharge;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class SimpleElectrostaticDirectMethod {
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(1);

    private val atoms : Rail[MMAtom];

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell, divided into separate Rails for each place
     */
    public def this(atoms : Rail[MMAtom]) {
		this.atoms = atoms;
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        var energy:Double = 0.0;

        for (i in 0..(atoms.size-1)) {
            val atomI = atoms(i);
            val ci = atomI.centre;
            val xi = ci.i;
            val yi = ci.j;
            val zi = ci.k;
            val qi = atomI.charge;
            var fix:Double = 0.0;
            var fiy:Double = 0.0;
            var fiz:Double = 0.0;
            for (j in 0..(atoms.size-1)) {
                if (i==j) continue;
	            val atomJ = atoms(j);
                val cj = atomJ.centre;

                val dx = cj.i-xi;
                val dy = cj.j-yi;
                val dz = cj.k-zi;

                val r2 = (dx*dx + dy*dy + dz*dz);
                val invR2 = 1.0 / r2;
                val qq = qi * atomJ.charge;
                val invR = Math.sqrt(invR2);

                val e = invR * qq;
                val forceScaling = e * invR2;
                energy += e;

                val fx = forceScaling * dx;
                val fy = forceScaling * dy;
                val fz = forceScaling * dz;

                fix += fx;
                fiy += fy;
                fiz += fz;
                //atomJ.force += Vector3d(-fx, -fy, -fz);
            }
            atomI.force = Vector3d(fix, fiy, fiz);
        }
       
        timer.stop(TIMER_INDEX_TOTAL);
        return 0.5 * energy;
    }

}
