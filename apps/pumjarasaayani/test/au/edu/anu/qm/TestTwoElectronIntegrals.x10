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
package au.edu.anu.qm;

import au.edu.anu.qm.TwoElectronIntegrals;

/**
 * Tests of base operations in the calculation of two-electron integrals
 * @author milthorpe
 */
public class TestTwoElectronIntegrals {

    public static def main(args : Rail[String]) {
	    val twoE = new TwoElectronIntegrals(3);
        val zeroM = twoE.computeZeroM(1, 2.5, 2.0, 0.5);
        for (i in 0..(zeroM.size-1)) {
            Console.OUT.println(i + " " + zeroM(i));
        }
    }
}

