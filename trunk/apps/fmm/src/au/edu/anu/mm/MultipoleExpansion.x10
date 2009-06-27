package au.edu.anu.mm;

import x10.lang.Complex;
import x10x.vector.Point3d;
import x10x.polar.Polar3d;

/**
 * This class calculates multipole expansions.
 * @author milthorpe
 */
public value MultipoleExpansion {

    /**
     * Calculate the multipole-like term O_{lm} (with m >= 0) for a point v.
     */
    public static def getOlm(q : double, v : Point3d, p : int) : Array[Complex]{rank==2} {
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
   * Calculate the Taylor-like term M_{lm} (with m >= 0) for a single point v.
   */
    public static def getMlm(v : Point3d, p : int) : Array[Complex]{rank==2} {
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
}


