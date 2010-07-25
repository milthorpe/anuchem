/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
package au.edu.anu.mm;

import harness.x10Test;

/**
 * Superclass for tests that require comparison of Complex and Double
 * @author milthorpe
 */
abstract class MathTest extends x10Test {

    public def nearlyEqual(a : Complex, b : Complex, maxRelativeError : Double, maxAbsoluteError : Double) : Boolean {
        if (a == b)
            return true;

        return nearlyEqual(a.re, b.re, maxRelativeError, maxAbsoluteError) 
            && nearlyEqual(a.im, b.im, maxRelativeError, maxAbsoluteError);
    }

    public def nearlyEqual(a: Double, b: Double, maxRelativeError : Double, maxAbsoluteError : Double) : Boolean {
        if (a == b)
            return true;

        if (Math.abs(a - b) < maxAbsoluteError)
            return true;

        var relativeError : Double;

        if (Math.abs(b) > Math.abs(a))
            relativeError = Math.abs((a - b) / b);
        else
            relativeError = Math.abs((a - b) / a);

        if (relativeError <= maxRelativeError)
            return true;
        else
            return false;
    }
}
