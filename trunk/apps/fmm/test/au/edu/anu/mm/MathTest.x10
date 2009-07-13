package au.edu.anu.mm;

import harness.x10Test;

/**
 * Superclass for tests that require comparison of Complex and Double
 * @author milthorpe
 */
abstract class MathTest extends x10Test {

    public def nearlyEqual(a : Complex, b : Complex, maxRelativeError : Double, maxAbsoluteError : Double) : Boolean {
        if (a.equals(b))
            return true;

        return nearlyEqual(a.real, b.real, maxRelativeError, maxAbsoluteError) 
            && nearlyEqual(a.imaginary, b.imaginary, maxRelativeError, maxAbsoluteError);
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
