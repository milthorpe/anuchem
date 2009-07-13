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

        val Mlm : LocalExpansion = LocalExpansion.getMlm(x, p);
        Console.OUT.println("local expansion:\n" + Mlm.toString());

        val target : LocalExpansion = new LocalExpansion(p);
        LocalExpansion.translateAndAddLocal(new Point3d(2.0, -3.0, 1.0), Mlm, target);
        Console.OUT.println("translated local:\n" + target.toString());

        val roundtrip : LocalExpansion = new LocalExpansion(p);
        LocalExpansion.translateAndAddLocal(new Point3d(-2.0, 3.0, -1.0), target, roundtrip);
        Console.OUT.println("translated local - roundtrip:\n" + roundtrip.toString());
		for (val(i,j) : Point in [0..p,-p..p]) {
            chk(nearlyEqual(roundtrip.terms(i,j), Mlm.terms(i,j), 1.0e-6, 1.0e-12));
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestLocalExpansion().execute();
    }

}
