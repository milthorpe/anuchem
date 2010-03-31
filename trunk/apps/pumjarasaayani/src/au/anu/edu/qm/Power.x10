/**
 * Power.x10
 *
 * Represents a gaussian power
 *
 * @author: V.Ganesh
 */
package au.anu.edu.qm;

public struct Power(l:Int, m:Int, n:Int) { 
    global val maxam:Int, minam:Int, totam:Int;

    public def this() { property(0,0,0); maxam=minam=totam=0; }

    public def this(l:Int, m:Int, n:Int) { 
        property(l, m, n);
        maxam = Math.max(l, Math.max(m, n));
        minam = Math.min(l, Math.min(m, n));
        totam = l+m+n;
    } 

    public def getL() = this.l;
    public def getM() = this.m;
    public def getN() = this.n;

    public def getTotalAngularMomentum() = totam;
    public def getMaximumAngularMomentum() = maxam;
    public def getMinimumAngularMomentum() = minam;
}

