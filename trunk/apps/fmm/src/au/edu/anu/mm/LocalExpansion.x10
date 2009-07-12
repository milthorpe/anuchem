package au.edu.anu.mm;

import x10.lang.Complex;
import x10x.vector.Tuple3d;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.polar.Polar3d;

/**
 * This class calculates local Taylor-type expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public value LocalExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }

    /**
     * Calculate the local Taylor-type expansion M_{lm} (with m >= 0) for a single point v.
     */
    public static def getMlm(v : Tuple3d, p : int) : LocalExpansion {
        val exp : LocalExpansion = new LocalExpansion(p);
        val v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double]{rank==2} = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 

        val rfac0 : Double = 1.0 / v_pole.r;
        val phifac0 : Complex = new Complex(Math.cos(v_pole.phi), Math.sin(v_pole.phi));
        var rfac : Double = rfac0;
        var il : int = 1;
        for (val (l): Point in [0..p]) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac.multiply(rfac * pplm(l,0) * ilm);
            for (val (m): Point in [1..l]) {
                ilm = ilm / (l+1-m);
                phifac = phifac.multiply(phifac0);
                exp.terms(l,m) = phifac.multiply(rfac * pplm(l,m) * ilm);
            }
            for (val (m): Point in [-l..-1]) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate().multiply((2*((-m+1)%2)-1 as Double));
            }
            rfac = rfac * rfac0;
        }

        return exp;
    }

    /** 
     * Translate a local expansion centred around the origin along a vector b,
     * and adds to a target expansion centred at b.
     * This corresponds to "Operator C", Equations 18-21 in White & Head-Gordon.
     * Note: this defines C^lm_jk(b) = O_j-l,k-m(b); however this is only defined
     * where abs(k-m) <= j-l, therefore we restrict k to [l-j+m..-l+j+m]
     * @param b the vector along which to translate the multipole
     * @param source the source local expansion, centred at the origin
     * @param target the target local expansion, centred at b
     */
    public static def translateAndAddLocal(b : Tuple3d,
                                         source : LocalExpansion,
                                         target : LocalExpansion) {
        val p : Int = source.terms.region.max(0);
        val shift : MultipoleExpansion = MultipoleExpansion.getOlm(b, p);
        for (val (l): Point in [0..p]) {
            for (val (j): Point in [l..p]) {
                for (val (m): Point in [-l..l]) {
                    for (val (k): Point in [l-j+m..-l+j+m]) {
                        val C_lmjk : Complex = shift.terms(j-l, k-m);
                        target.terms(l,m) = target.terms(l,m).add(C_lmjk.multiply(source.terms(j,k)));
                    }
                }
            }
        }
    }

}


