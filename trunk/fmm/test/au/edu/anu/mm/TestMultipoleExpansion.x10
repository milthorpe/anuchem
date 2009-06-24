package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;
import au.edu.anu.mm.MultipoleExpansion;

/**
 * Test calculation of associated Legendre polynomials.
 * @author milthorpe
 */
class TestMultipoleExpansion extends x10Test {
    public def run(): boolean {
        val x : Point3d = new Point3d(1.0, 2.0, -1.0);
        val Mlm : Array[Complex]{rank==2} = MultipoleExpansion.getMlm(x, 3);
		for (val(i) : Point in [0..3]) {
            for (val(j) : Point in [-3..3]) {
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
