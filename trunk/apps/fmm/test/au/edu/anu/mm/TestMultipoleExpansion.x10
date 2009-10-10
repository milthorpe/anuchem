package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;

/**
 * Test multipole expansions
 * @author milthorpe
 */
class TestMultipoleExpansion extends MathTest {
    public def run(): boolean {
        val p : Int = 3; // multipole level

        val x : Point3d = new Point3d(1.0, 2.0, -1.0);
        val Olm : MultipoleExpansion = MultipoleExpansion.getOlm(1.5, x, p);
        Console.OUT.println("multipole expansion:\n" + Olm.toString());

        val target : MultipoleExpansion = new MultipoleExpansion(p);
        val translation = MultipoleExpansion.getOlm(new Point3d(2.0, -3.0, 1.0), p);
        target.translateAndAddMultipole(translation, Olm);
        Console.OUT.println("translated multipole:\n" + target.toString());

        val roundtrip : MultipoleExpansion = new MultipoleExpansion(p);
        val reverseTranslation = MultipoleExpansion.getOlm(new Point3d(-2.0, 3.0, -1.0), p);
        roundtrip.translateAndAddMultipole(reverseTranslation, target);
        Console.OUT.println("translated multipole - roundtrip:\n" + roundtrip.toString());
		for ((i): Point in [0..p]) {
            for ((j): Point in [-i..i]) {
                chk(nearlyEqual(roundtrip.terms(i,j), Olm.terms(i,j), 1.0e-6, 1.0e-12)); 
		    }
        }

        val localExp : LocalExpansion = new LocalExpansion(p);
        val transformation = LocalExpansion.getMlm(new Point3d(2.0, -3.0, 1.0), p);
        localExp.transformAndAddToLocal(transformation, Olm);
        Console.OUT.println("transformed multipole:\n" + localExp.toString());

        return true;
    }

    public static def main(Rail[String]) {
        new TestMultipoleExpansion().execute();
    }

}
