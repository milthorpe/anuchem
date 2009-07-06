package au.edu.anu.mm;

import x10.lang.Complex;
import x10x.vector.Tuple3d;
import x10x.vector.Point3d;
import x10x.vector.Vector3d;
import x10x.polar.Polar3d;

/**
 * This class calculates multipole expansions, using the algorithms given in
 * <a href="info:doi/10.1063/1.468354">
 *      Derivation and efficient implementation of the fast multipole method
 * </a>, White and Head-Gordon, 1994.
 * @author milthorpe
 */
public value MultipoleExpansion extends Expansion {

    public def this(p : Int) {
        super(p);
    }

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(q : Double, v : Tuple3d, p : int) : MultipoleExpansion {
        val exp : MultipoleExpansion = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double]{rank==2} = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 : Complex = new Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Int = 1;
        for (val (l): Point in [0..p]) {
            il = il * Math.max(l,1);
            var ilm : Int = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac.divide(ilm).multiply(q * rfac * pplm(l,0)); 
            for (val (m): Point in [1..l]) {
                ilm = ilm*(l+m);
                phifac = phifac.multiply(phifac0);
                exp.terms(l,m) = phifac.divide(ilm).multiply(q * rfac * pplm(l,m));
            }
            for (val (m): Point in [-l..-1]) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate().multiply((2*((-m+1)%2)-1 as Double));
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /**
     * Calculate the chargeless multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(v : Tuple3d, p : int) : MultipoleExpansion {
        val exp : MultipoleExpansion = new MultipoleExpansion(p);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[Double]{rank==2} = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 : Complex = new Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : Double = 1.0;
        var il : Int = 1;
        for (val (l): Point in [0..p]) {
            il = il * Math.max(l,1);
            var ilm : Int = il;
            var phifac : Complex = Complex.ONE;
            exp.terms(l,0) = phifac.divide(ilm).multiply(rfac * pplm(l,0)); 
            for (val (m): Point in [1..l]) {
                ilm = ilm*(l+m);
                phifac = phifac.multiply(phifac0);
                exp.terms(l,m) = phifac.divide(ilm).multiply(rfac * pplm(l,m));
            }
            for (val (m): Point in [-l..-1]) {
                exp.terms(l,m) = exp.terms(l,-m).conjugate().multiply((2*((-m+1)%2)-1 as Double));
            }
            rfac = rfac * v_pole.r;
        }

        return exp;
    }

    /** 
     * Translate a multipole expansion centred around the origin along a vector -b,
     * and adds to a target expansion centred at -b.
     * This corresponds to "Operator A", Equations 10-11 in White & Head-Gordon.
     * Note: this defines A^lm_jk(b) = O_l-j,m-k(b); however this is only defined
     * where abs(m-k) <= l-j, therefore we restrict m to [j-l+k..-j+l+k]
     * @param b the vector along which to translate the multipole
     * @param source the source multipole expansion, centred at the origin
     * @param target the target multipole expansion, centred at -b
     */
    public static def translateAndAddMultipole(b : Tuple3d,
                                         source : MultipoleExpansion,
                                         target : MultipoleExpansion) {
        val p : Int = source.terms.region.max(0);
        val shift : MultipoleExpansion = MultipoleExpansion.getOlm(b, p);
        for (val (j): Point in [0..p]) {
            for (val (l): Point in [j..p]) {
                for (val (k): Point in [-j..j]) {
                    for (val (m): Point in [j-l+k..-j+l+k]) {
                        val A_lmjk : Complex = shift.terms(l-j, m-k);
                        target.terms(l,m) = target.terms(l,m).add(A_lmjk.multiply(source.terms(j,k)));
                    }
                }
            }
        }
    }

    /**
     * Perform local expansion translation along a vector V.
     * This is used in FMM to translate from parent to child box.
     * // TODO use correct array regions (not shoehorned into 1D array)
     * @param v the vector from the centre of parent box to centre of child box
     * @param parentExpansion the parent local expansion
     * @param childExpansion the child local expansion to which the parent will be added
     */
    public static def translateExpansion(v : Tuple3d, 
                            parentExpansion: Array[Complex]{rank==1},
                            var childExpansion : Array[Complex]{rank==1}) {

        val numTerms : Int = (Math.sqrt(parentExpansion.region.max(0)) as Int) - 1;
        val cylo : Array[Complex]{rank==2} = MultipoleExpansion.getOlm(1.0, v, numTerms).terms;
        var inm : Int = 0;
        for (val (jt): Point in [0..numTerms]) {
            for (val (kt) : Point in [-jt..jt]) {
                inm++;
                for (val(lt): Point in [jt..numTerms]) {
                    for (val (mt): Point in [kt-(lt-jt)..kt+(lt-jt)]) {
                        //Console.OUT.println("jt = " + jt + " kt = " + kt + " lt = " + lt + " mt = " + mt);  
                        val ido = (lt*lt)+lt+mt+1;
                       // Console.OUT.println("inm = " + inm + " ido = " + ido);
                        childExpansion(inm) = childExpansion(inm).add(parentExpansion(ido).multiply(cylo(lt-jt,mt-kt)));
                    }
                }
            }
        }
    }
}


