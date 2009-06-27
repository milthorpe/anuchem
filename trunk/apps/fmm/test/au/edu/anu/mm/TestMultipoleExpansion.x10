package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;
import au.edu.anu.mm.MultipoleExpansion;

/**
 * Test multipole expansion including Taylor-like terms.
 * @author milthorpe
 */
class TestMultipoleExpansion extends x10Test {
    public def run(): boolean {
        val p : int = 2; // multipole level

        val x : Point3d = new Point3d(1.0, 2.0, -1.0);
        val Olm : Array[Complex]{rank==2} = MultipoleExpansion.getOlm(1.5, x, p);
        Console.OUT.println("multipole expansion");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(Olm(i,j) + " ");
            }
            Console.OUT.println();
		}

        val Mlm : Array[Complex]{rank==2} = MultipoleExpansion.getMlm(x, p);
        Console.OUT.println("local expansion");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(Mlm(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestMultipoleExpansion().execute();
    }

}
