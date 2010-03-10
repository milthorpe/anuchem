/**
 * Power.x10
 *
 * Represents a gaussian power
 *
 * @author: V.Ganesh
 */
package au.anu.edu.qm;

public struct Power(l:Int, m:Int, n:Int) { 
    public def this() { property(0,0,0); }

    public def this(l:Int, m:Int, n:Int) { 
        property(l, m, n);
    } 

    public def getL() : Int = this.l;
    public def getM() : Int = this.m;
    public def getN() : Int = this.n;

    public def getTotalAngularMomentum() = l+m+n;
    public def getMaximumAngularMomentum() = Math.max(l, Math.max(m, n));
    public def getMinimumAngularMomentum() = Math.min(l, Math.min(m, n));
}

