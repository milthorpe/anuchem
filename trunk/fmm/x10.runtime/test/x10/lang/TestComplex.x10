package x10.lang;

import harness.x10Test;

/**
 * @author milthorpe
 */
class TestComplex extends x10Test {
    public def run(): boolean {
        val a : Complex = new Complex(2.0, 2.0);
		chk (a.negate().negate() == a);
		chk (a.abs() == Math.sqrt(8.0));

		val b : Complex = a.conjugate();
		chk (b.conjugate() == a);

		val c : Complex = new Complex(1.0, 3.0);
		chk (a.add(c).subtract(c) == a, "a + c - c = a");
		chk (a.multiply(c).divide(c) == a, "a * c / c = a");

		val d : Complex = new Complex(4.0, -1.0);
		chk (a.add(d).subtract(d) == a, "a + d - d = a");
		chk (a.multiply(d).divide(d) == a, "a * d / d = a");

		val e : Complex = a.divide(Complex.ZERO);
		chk (e.isNaN());

        return true;
    }

    public static def main(Rail[String]) {
        new TestComplex().execute();
    }

}
