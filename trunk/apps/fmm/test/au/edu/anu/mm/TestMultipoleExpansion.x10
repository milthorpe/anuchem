package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;

/**
 * Test multipole expansions
 * @author milthorpe
 */
class TestMultipoleExpansion extends x10Test {
    public def run(): boolean {
        val p : Int = 1; // multipole level

        val x : Point3d = new Point3d(1.0, 2.0, -1.0);
        val Olm : MultipoleExpansion = MultipoleExpansion.getOlm(1.5, x, p);
        Console.OUT.println("multipole expansion");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(Olm.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        val target : MultipoleExpansion = new MultipoleExpansion(p);
        MultipoleExpansion.translateAndAddMultipole(new Point3d(2.0, -3.0, 1.0), Olm, target);
        Console.OUT.println("translated multipole");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(target.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        val roundtrip : MultipoleExpansion = new MultipoleExpansion(p);
        MultipoleExpansion.translateAndAddMultipole(new Point3d(-2.0, 3.0, -1.0), target, roundtrip);
        Console.OUT.println("translated - roundtrip");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(roundtrip.terms(i,j) + " ");
                //chk(nearlyEqual(roundtrip.terms(i,j), Olm.terms(i,j), 1.0e-6));
            }
            Console.OUT.println();
		}
/*
        val localExp : LocalExpansion = new LocalExpansion(p);
        MultipoleExpansion.transformAndAddToLocal(new Point3d(2.0, -3.0, 1.0), Olm, localExp);
        Console.OUT.println("transformed multipole");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(localExp.terms(i,j) + " ");
            }
            Console.OUT.println();
		}
*/
        return true;
    }


    public def nearlyEqual(a : Complex, b : Complex, maxRelativeError : Double) : Boolean {
        if (a.equals(b))
            return true;

        return nearlyEqual(a.real, b.real, maxRelativeError) && nearlyEqual(a.imaginary, b.imaginary, maxRelativeError);
    }

    public def  nearlyEqual(a: Double, b: Double, maxRelativeError : Double) : Boolean {
        if (a == b)
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

    public static def main(Rail[String]) {
        new TestMultipoleExpansion().execute();
    }

}
