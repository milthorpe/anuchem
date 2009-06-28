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
public value MultipoleExpansion {

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(q : double, v : Tuple3d, p : int) : Array[Complex]{rank==2} {
        val Olm : Array[Complex]{rank==2} = Array.make[Complex]([0..p, -p..p]->here, (val (i,j): Point)=> Complex.ZERO);
        var v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[double]{rank==2} = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 
        
        val phifac0 : Complex = new Complex(Math.cos(-v_pole.phi), Math.sin(-v_pole.phi));

        var rfac : double = 1.0;
        var il : int = 1;
        for (val (l): Point in [0..p]) {
            il = il * Math.max(l,1);
            var ilm : int = il;
            var phifac : Complex = Complex.ONE;
            Olm(l,0) = phifac.divide(ilm).multiply(q * rfac * pplm(l,0)); 
            for (val (m): Point in [1..l]) {
                ilm = ilm*(l+m);
                phifac = phifac.multiply(phifac0);
                Olm(l,m) = phifac.divide(ilm).multiply(q * rfac * pplm(l,m));
            }
            for (val (m): Point in [-l..-1]) {
                Olm(l,m) = Olm(l,-m).conjugate().multiply((2*((-m+1)%2)-1 as double));
            }
            rfac = rfac * v_pole.r;
        }

        return Olm;
    }

    /**
     * Calculate the local Taylor-type expansion M_{lm} (with m >= 0) for a single point v.
     */
    public static def getMlm(v : Tuple3d, p : int) : Array[Complex]{rank==2} {
        val Mlm : Array[Complex]{rank==2} = Array.make[Complex]([0..p, -p..p]->here, (val (i,j): Point)=> Complex.ZERO);
        val v_pole : Polar3d = Polar3d.getPolar3d(v);
        val pplm : Array[double]{rank==2} = AssociatedLegendrePolynomial.getPlm(Math.cos(v_pole.theta), p); 

        val rfac0 : double = 1.0 / v_pole.r;
        val phifac0 : Complex = new Complex(Math.cos(v_pole.phi), Math.sin(v_pole.phi));
        var rfac : double = rfac0;
        var il : int = 1;
        for (val (l): Point in [0..p]) {
            il = il * Math.max(l,1);
            var ilm : double = il;
            var phifac : Complex = Complex.ONE;
            Mlm(l,0) = phifac.multiply(rfac * pplm(l,0) * ilm);
            for (val (m): Point in [1..l]) {
              ilm = ilm / (l+1-m);
              phifac = phifac.multiply(phifac0);
              Mlm(l,m) = phifac.multiply(rfac * pplm(l,m) * ilm);
            }
            for (val (m): Point in [-l..-1]) {
              Mlm(l,m) = Mlm(l,-m).conjugate().multiply((2*((-m+1)%2)-1 as double));
            }
            rfac = rfac * rfac0;
        }

        return Mlm;
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

        val numTerms : int = parentExpansion.region.max(0);
        val cylo : Array[Complex]{rank==2} = getOlm(1.0, v, numTerms);
        var inm : int = 0;
        for (val (jt): Point in [0..numTerms]) {
            for (val (kt) : Point in [jt..jt]) {
                inm++;
                for (val(lt): Point in [jt..numTerms]) {
                    for (val (mt): Point in [kt-(lt-jt)..kt+(lt-jt)]) {     
                        val ido = (lt*lt)+lt+mt+1;
                        childExpansion(inm) = childExpansion(inm).add(parentExpansion(ido).multiply(cylo(lt-jt,mt-kt)));
                    }
                }
            }
        }
    }
}


