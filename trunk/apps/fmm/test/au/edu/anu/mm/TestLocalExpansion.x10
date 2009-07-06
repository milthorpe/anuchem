package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;

/**
 * Test local Taylor-type expansions.
 * @author milthorpe
 */
class TestLocalExpansion extends x10Test {
    public def run(): boolean {
        val p : int = 3; // multipole level
        val x : Point3d = new Point3d(1.0, 2.0, -1.0);

        val Mlm : LocalExpansion = LocalExpansion.getMlm(x, p);
        Console.OUT.println("local expansion");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(Mlm.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        val target : LocalExpansion = new LocalExpansion(p);
        LocalExpansion.translateAndAddLocal(new Point3d(2.0, -3.0, 1.0), Mlm, target);
        Console.OUT.println("translated local");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(target.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestLocalExpansion().execute();
    }

}
