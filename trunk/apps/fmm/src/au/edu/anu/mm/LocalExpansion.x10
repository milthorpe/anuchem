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

    public def this(p : Int, data : ValRail[Complex]) {
        super(p, data);
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
        for ((l) in 0..p) {
            il = il * Math.max(l,1);
            var ilm : Double = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac * (rfac * pplm(l,0) * ilm);
            for ((m) in 1..l) {
                ilm = ilm / (l+1-m);
                phifac = phifac * phifac0;
                exp.terms(l,m) = phifac * (rfac * pplm(l,m) * ilm);
            }
            for ((m) in -l..-1) {
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
    public safe def translateAndAddLocal(shift : MultipoleExpansion!,
                                         source : LocalExpansion!) {
        val p = terms.region.max(0);
        for ((l,m) in terms.dist) {
            for ((j) in l..p) {
                for ((k) in (l-j+m)..(-l+j+m)) {
                    val C_lmjk = shift.terms(j-l, k-m);
                    val O_jk = source.terms(j,k);
                    this.terms(l,m) = this.terms(l,m) + C_lmjk * O_jk;
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
    public safe def transformAndAddToLocal(transform : LocalExpansion!,
                                         source : MultipoleExpansion!) {
        val p : Int = terms.region.max(0);
        for ((j,k) in terms.region) {
            val O_jk = source.terms(j,k);
            for ((l) in 0..p-j) {
                for ((m) in -l..l) {
                    if (Math.abs(k+m) <= (j+l)) {
                        val B_lmjk = transform.terms(j+l, k+m);
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
    public safe def getPotential(q : Double,
                                v : Tuple3d) : Double {
        val numTerms = terms.region.max(0);
        val transform = MultipoleExpansion.getOlm(q, v, numTerms);
        var potential : Double = 0.0;
        // TODO use lift/reduction?
        for ((i,j) in terms.region) {
            potential += (terms(i,j) * transform.terms(i,j)).re;
        }
        return potential;
    }

    /**
     * TODO this should not be necessary once XTENLANG-787 is resolved
     * @return a local copy of a local expansion from another place
     */
    static def getLocalCopy(p : Int, source : LocalExpansion) : LocalExpansion! {
        val data = at (source) {getData(p, source)};
        return new LocalExpansion(p, data);
    }
}


