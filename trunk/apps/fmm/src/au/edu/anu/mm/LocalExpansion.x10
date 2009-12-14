package au.edu.anu.mm;

import x10x.vector.Tuple3d;
import x10x.polar.Polar3d;

/**
 * This class calculates local Taylor-type expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public class LocalExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }

    /**
     * Calculate the local Taylor-type expansion M_{lm} (with m >= 0) for a single point v.
     */
    public static safe def getMlm(v : Tuple3d, p : int) : LocalExpansion! {
        val exp = new LocalExpansion(p);
        val v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double](2) = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 

        val rfac0 : Double = 1.0 / v_pole.r;
        val phifac0 = Complex(Math.cos(v_pole.phi), Math.sin(v_pole.phi));
        var rfac : Double = rfac0;
        var il : Double = 1.0;
        for (var l : Int = 0; l<=p; l++) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac * (rfac * pplm(l,0) * ilm);
            for (var m : Int = 1; m<=l; m++) {
                ilm = ilm / (l+1-m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac * (rfac * pplm(l,m) * ilm);
            }
            for (var m : Int = -l; m<=-1; m++) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate() * ((2*((-m+1)%2)-1 as Double));
            }
            rfac = rfac * rfac0;
        }

        return exp;
    }

    /** 
     * Translate a local expansion centred around the origin along a vector b,
     * and adds to this expansion centred at b, where shift is the expansion
     * of the vector b.
     * This corresponds to "Operator C", Equations 18-21 in White & Head-Gordon.
     * Note: this defines C^lm_jk(b) = O_j-l,k-m(b); however this is only defined
     * where abs(k-m) <= j-l, therefore we restrict k to [l-j+m..-l+j+m]
     * @param shift the multipole expansion of the translation
     * @param source the source local expansion, centred at the origin
     */
    public global safe def translateAndAddLocal(shift : MultipoleExpansion!,
                                         source : LocalExpansion!) {
        val p : Int = source.terms.region.max(0);
        for (val (l,m): Point in this.terms) {
            for (var j : Int = l; j<=p; j++) { // TODO XTENLANG-504
                for (var k : Int = l-j+m; k<=-l+j+m; k++) { // TODO XTENLANG-504
                    val C_lmjk = shift.terms(j-l, k-m);
                    this.terms(l,m) = this.terms(l,m) + C_lmjk * source.terms(j,k);
                }
            }
        }
    }

   /** 
     * Transform a multipole expansion centred around the origin into a
     * Taylor expansion centred about b, and adds to this expansion.
     * This corresponds to "Operator B", Equations 13-15 in White & Head-Gordon.
     * This operator is inexact due to truncation of the series at <em>p</em> poles.
     * Note: this defines B^lm_jk(b) = M_j+l,k+m(b), therefore restrict l to [0..p-j]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     */
    public global def transformAndAddToLocal(transform : LocalExpansion,
                                         source : MultipoleExpansion) {
        val p : Int = source.terms.region.max(0);
        for (val (j,k): Point in source.terms) {
            val O_jk = source.terms(j,k);
            for (var l : Int = 0; l <= p-j; l++) { // TODO XTENLANG-504
                for (var m : Int = -l; m<=l; m++) { // TODO XTENLANG-504
                    if (Math.abs(k+m) <= (j+l)) {
                        // TODO calculating the indices "on the fly" in the body 
                        // of the at statement results in a Seg Fault... why?
                        val jPlusL : Int = (j+l);
                        val kPlusM : Int = (k+m);
                        val B_lmjk = transform.terms(jPlusL, kPlusM);
                        //Console.OUT.println("source.terms.dist(" + j + "," + k + ") = " + source.terms.dist(j,k));
                        this.terms(l,m) = this.terms(l,m) + B_lmjk * O_jk;
                    }
                }
            }
        }
    }

    /**
     * Transforms this local expansion about the origin to the potential
     * acting on <code>q</code> at point <code>v</code>.
     * @param q the charge at point v
     * @param v the location of charge q
     */
    public global safe def getPotential(q : Double,
                                v : Tuple3d) : Double {
        val numTerms = terms.region.max(0);
        val transform = MultipoleExpansion.getOlm(q, v, numTerms);
        var potential : Double = 0.0;
        // TODO use lift/reduction?
        for (p in terms.region) {
            potential += (terms(p) * transform.terms(p)).re;
        }
        return potential;
    }
}


