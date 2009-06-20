package au.edu.anu.mm;

import harness.x10Test;
import au.edu.anu.mm.AssociatedLegendrePolynomial;

/**
 * Test calculation of associated Legendre polynomials.
 * @author milthorpe
 */
class TestAssociatedLegendrePolynomial extends x10Test {
    public def run(): boolean {
        val Plm : Array[double]{rank==2} = AssociatedLegendrePolynomial.getPlm(0.5, 3);
		for (val(i) : Point in [0..3]) {
            for (val(j) : Point in [0..i]) {
			    Console.OUT.print(Plm(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestAssociatedLegendrePolynomial().execute();
    }

}
