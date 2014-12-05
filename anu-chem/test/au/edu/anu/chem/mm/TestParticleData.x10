/*
 *  This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 *  (C) Copyright IBM Corporation 2014.
 */
package au.edu.anu.chem.mm;

import x10.util.Random;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

import harness.x10Test;

public class TestParticleData extends x10Test {
    def testSort():Boolean {
        val rand = new Random();
        val particleData = new ParticleData();
        for (i in 0..9) {
            particleData.addAtom(i, 0n, new Point3d(rand.nextDouble(), rand.nextDouble(), rand.nextDouble()), Vector3d.NULL);
        }
        val unsortedCenters = particleData.x.toRail();
/*
        for (i in 0..(particleData.numAtoms()-1)) {
            Console.OUT.println(particleData.globalIndex(i) + " " + particleData.x(i));
        }
*/
        particleData.sortAtoms(0, particleData.numAtoms()-1, (a:Long,b:Long)=>particleData.x(a).i.compareTo(particleData.x(b).i));

        for (i in 1..(particleData.numAtoms()-1)) {
            // check particles are in order of particle center i coordinate
            chk(particleData.x(i).i >= particleData.x(i-1).i);
            // and that sort preserved mapping between fields
            chk(unsortedCenters(particleData.globalIndex(i)) == particleData.x(i));
        }

        return true;
    }

    public def run(): boolean {
        var result:Boolean = true;
        result &= testSort();

        return result;
    }

    public static def main(args:Rail[String]) {
        new TestParticleData().execute();
    }
}
