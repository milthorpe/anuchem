package au.edu.anu.mm;

/**
 * This class calculates associated Legendre polynomials.
 * @author milthorpe
 */
public value AssociatedLegendrePolynomial {

	/**
	 * Calculate associated Legendre polynomials P_{lm}(x) up to l=p (m.ge.0)
	 */
	public static def getPlm(x: double, p : int) : Array[double]{rank==2} {
        if (Math.abs(x) > 1.0) {
            throw new IllegalArgumentException("abs(x) > 1: Associated Legendre functions are only defined on [-1, 1].");
        }

		val Plm : Array[double]{rank==2} = Array.make[double](Region.makeLowerTriangular(p+1)->here);
		Plm(0,0) = 1.0;
		val somx2 : double = Math.sqrt((1.0 - x) * (1.0 + x));
		var fact : double = 1.0;
		for (val (i): Point in [1..p]) {
			Plm(i,i) = -Plm(i-1,i-1) * fact * somx2;
			Plm(i,i-1) = x * fact * Plm(i-1,i-1);
			fact += 2.0;
		}
		for (val (m): Point in [0..p-2]) {
			for (val (l): Point in [m+2..p]) {
				Plm(l,m) = (x * (2.0 * l - 1.0) 
							* Plm(l-1,m) - (l + m - 1.0) 
							* Plm(l-2,m)) / (l - m);
			}
		}
		return Plm;
	}

}