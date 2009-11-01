package au.edu.anu.mm;

import x10x.vector.Tuple3d;
import x10x.polar.Polar3d;

/**
 * This class calculates multipole expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public class MultipoleExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     * TODO x10doc
     */
    public static def getOlm(q : Double, v : Tuple3d, p : Int) : MultipoleExpansion! {
        val exp = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double](2) = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Double = 1.0;
        for (var l : Int = 0; l<=p; l++) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac / ilm * (q * rfac * pplm(l,0)); 
            for (var m : Int = 1; m<=l; m++) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac / ilm * (q * rfac * pplm(l,m));
            }
            for (var m : Int = -l; m<=-1; m++) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate() * (2*((-m+1)%2)-1);
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /**
     * Calculate the chargeless multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(v : Tuple3d, p : Int) : MultipoleExpansion! {
        val exp = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double](2) = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 = Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Double = 1.0;
        for (var l : Int = 0; l<=p; l++) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac / ilm * (rfac * pplm(l,0)); 
            for (var m : Int = 1; m<=l; m++) {
                ilm = ilm*(l+m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac / ilm * (rfac * pplm(l,m));
            }
            for (var m : Int = -l; m<=-1; m++) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate() * (2*((-m+1)%2)-1);
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /** 
     * Translate a multipole expansion centred around the origin along a vector -b,
     * and adds to this expansion centred at -b.
     * This corresponds to "Operator A", Equations 10-11 in White & Head-Gordon.
     * Note: this defines A^lm_jk(b) = O_l-j,m-k(b); however this is only defined
     * where abs(m-k) <= l-j, therefore we restrict m to [j-l+k..-j+l+k]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     */
    public def translateAndAddMultipole(shift : MultipoleExpansion!,
                                         source : MultipoleExpansion!) {
        val p : Int = source.terms.region.max(0);
        for (val (j,k): Point in source.terms) {
            val O_jk = source.terms(j,k);
            for (var l : Int = j; l<=p; l++) { // TODO XTENLANG-504
                for (var m : Int = -l; m<=l; m++) { // TODO XTENLANG-504
                    if (Math.abs(m-k) <= (l-j)) {
                        val A_lmjk = shift.terms(l-j, m-k);
                        this.terms(l,m) = this.terms(l,m) + A_lmjk * O_jk;
                    }
                }
            }
        }
    }

    /** 
     * Translate a multipole expansion centred around the origin along a vector -b,
     * and adds to this expansion centred at -b.
     * This corresponds to "Operator A", Equations 10-11 in White & Head-Gordon.
     * Note: this defines A^lm_jk(b) = O_l-j,m-k(b); however this is only defined
     * where abs(m-k) <= l-j, therefore we restrict m to [j-l+k..-j+l+k]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     */
    public def translateAndAddMultipoleDist(shift : MultipoleExpansion!,
                                         source : MultipoleExpansion) {
        val p : Int = this.terms.region.max(0);
        for (val (j,k): Point in this.terms) {
            val O_jk = at (source) {source.terms(j,k)};
            for (var l : Int = j; l<=p; l++) { // TODO XTENLANG-504
                for (var m : Int = -l; m<=l; m++) { // TODO XTENLANG-504
                    if (Math.abs(m-k) <= (l-j)) {
                        val A_lmjk = shift.terms(l-j, m-k);
                        this.terms(l,m) = this.terms(l,m) + A_lmjk * O_jk;
                    }
                }
            }
        }
    }
}


