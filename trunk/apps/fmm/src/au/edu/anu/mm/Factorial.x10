/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Andrew Haigh 2011.
 */

package au.edu.anu.mm;

/**
 * This class calculates factorials (up to a given limit) and stores them to avoid having to calculate in time-crucial operations
 * @author haigh
 */
public class Factorial { 
	public static val fact : Array[Double](1) = new Array[Double](200);

        /**
         * This function must be called first by every class that wants to
         * do a translation using rotations or generate Wigner matrices pre-multiplied (i.e. with WignerRotationMatrix.getExpandedCollection)
         */
	public static def calcFactorial(numTerms : int) { 
		fact(0) = 1.0;
		for ([i] in 1..(2*numTerms)) fact(i) = i * fact(i-1);
	}

	public static def getFactorial(i : int) = fact(i); 
}
