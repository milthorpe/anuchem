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
	public static val factorial : Array[Double](1) = Factorial.calcFact();

    private static final def calcFact() { 
        val fact = new Array[Double](100);
		fact(0) = 1.0;
		for ([i] in 1..99) fact(i) = i * fact(i-1);
        return fact;
	}

	public static def getFactorial(i : int) = factorial(i); 
}
