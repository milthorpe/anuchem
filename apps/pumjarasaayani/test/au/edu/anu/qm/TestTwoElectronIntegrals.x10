/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Josh Milthorpe 2011-2013.
 */
package au.edu.anu.qm;

import au.edu.anu.qm.TwoElectronIntegrals;

/**
 * Tests of base operations in the calculation of two-electron integrals
 * @author milthorpe
 */
public class TestTwoElectronIntegrals {

    public static def main(args : Rail[String]) {
        val norm = [1.0, 2.0, 1.0];
	    val twoE = new TwoElectronIntegrals(3n, norm, 0.1, 0.00001);
        val zeroM = twoE.computeZeroM(1n, 2.5, 2.0, 0.5);
        for (i in 0..(zeroM.size-1)) {
            Console.OUT.println(i + " " + zeroM(i));
        }
    }
}

