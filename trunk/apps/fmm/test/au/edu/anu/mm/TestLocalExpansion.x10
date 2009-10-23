package au.edu.anu.mm;

import x10x.vector.Point3d;

/**
 * Test local Taylor-type expansions.
 * @author milthorpe
 */
class TestLocalExpansion extends MathTest {
    public def run(): boolean {
        val p : int = 3; // multipole level
        val x : Point3d = new Point3d(1.0, 2.0, -1.0);

        val Mlm = LocalExpansion.getMlm(x, p);
        Console.OUT.println("local expansion:\n" + Mlm.toString());

        val target = new LocalExpansion(p);
        val translation = MultipoleExpansion.getOlm(new Point3d(2.0, -3.0, 1.0), p);
        target.translateAndAddLocal(translation, Mlm);
        Console.OUT.println("translated local:\n" + target.toString());

        val roundtrip = new LocalExpansion(p);
        val reverseTranslation = MultipoleExpansion.getOlm(new Point3d(-2.0, 3.0, -1.0), p);
        roundtrip.translateAndAddLocal(reverseTranslation, target);
        Console.OUT.println("translated local - roundtrip:\n" + roundtrip.toString());
		for ((i): Point in [0..p]) {
            for ((j): Point in [-i..i]) {
                chk(nearlyEqual(roundtrip.terms(i,j), Mlm.terms(i,j), 1.0e-6, 1.0e-12));
		    }
        }

        return true;
    }

    public static def main(Rail[String]) {
        new TestLocalExpansion().execute();
    }

}
