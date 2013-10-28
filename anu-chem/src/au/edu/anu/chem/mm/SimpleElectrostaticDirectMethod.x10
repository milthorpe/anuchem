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

import au.edu.anu.util.Timer;

/**
 * This class calculates electrostatic interactions between
 * particles directly.  This is an O(N^2) calculation, intended
 * for comparison with other methods e.g. FMM, SPME.
 */
public class SimpleElectrostaticDirectMethod {
    public static val SIZE = 80.0;
    // TODO enum - XTENLANG-1118
    public static val TIMER_INDEX_TOTAL = 0;
    /** A multi-timer for the several segments of a single getEnergy invocation, indexed by the constants above. */
    public val timer = new Timer(1);

    public static class Atom {
        public var x:Double;
        public var y:Double;
        public var z:Double;
        public var charge:Double;
        public var fx:Double;
        public var fy:Double;
        public var fz:Double;
    }

    private val atoms : Rail[Atom];

    /**
     * Creates a new electrostatic direct method.
     * @param atoms the atoms in the unit cell, divided into separate Rails for each place
     */
    public def this(atoms : Rail[Atom]) {
		this.atoms = atoms;
    }
	
    public def getEnergy() : Double {
        timer.start(TIMER_INDEX_TOTAL);

        var energy:Double = 0.0;

        for (i in 0..(atoms.size-1)) {
            val atomI = atoms(i);
            val xi = atomI.x;
            val yi = atomI.y;
            val zi = atomI.z;
            val qi = atomI.charge;
            var fix:Double = 0.0;
            var fiy:Double = 0.0;
            var fiz:Double = 0.0;
            for (j in 0..(atoms.size-1)) {
                if (i==j) continue;
	            val atomJ = atoms(j);

                val dx = atomJ.x-xi;
                val dy = atomJ.y-yi;
                val dz = atomJ.z-zi;

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
            atomI.fx += fix;
            atomI.fy += fiy;
            atomI.fz += fiz;
        }
       
        timer.stop(TIMER_INDEX_TOTAL);
        return 0.5 * energy;
    }

    public static def generateAtoms(numAtoms:Long):Rail[Atom] {
	    val atoms = new Rail[Atom](numAtoms, (Long) => new Atom());
        val gridSize = (Math.ceil(Math.cbrt(numAtoms))) as Int;
        var gridPoint:Long = 0; // running total of assigned grid points
	    for (i in 0..(numAtoms-1)) {
		    val gridX = gridPoint / (gridSize * gridSize);
            val gridY = (gridPoint - (gridX * gridSize * gridSize)) / gridSize;
            val gridZ = gridPoint - (gridX * gridSize * gridSize) - (gridY * gridSize);
            val x = (gridX + 0.5) * (SIZE / gridSize);
            val y = (gridY + 0.5) * (SIZE / gridSize);
            val z = (gridZ + 0.5) * (SIZE / gridSize);
            atoms(i).x = x;
		    atoms(i).y = y;
		    atoms(i).z = z;
		    atoms(i).charge = gridPoint%2==0?1.0:-1.0;
            atoms(i).fx = 0.0;
            atoms(i).fy = 0.0;
            atoms(i).fz = 0.0;
            gridPoint++;
	    }
        return atoms;
    }

}
