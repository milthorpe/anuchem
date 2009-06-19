package au.edu.anu.mm;

import harness.x10Test;
import au.edu.anu.mm.AssociatedLegendrePolynomial;

/**
 * Test calculation of associated Legendre polynomials.
 * @author milthorpe
 */
class TestAssociatedLegendrePolynomial extends x10Test {
    public def run(): boolean {
        val Plm : Array[double]{rank==2} = AssociatedLegendrePolynomial.getPlm(1.5, 2);
		for (val(i,j) : Point in Plm.region) {
			Console.OUT.println(Plm(i,j));
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestAssociatedLegendrePolynomial().execute();
    }

}
