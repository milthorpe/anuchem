/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010-2012.
 */
package au.edu.anu.qm;

import x10.compiler.NonEscaping;
import x10.util.ArrayList;
import x10x.vector.Point3d;

/**
 * Represents a list of ShellPair
 *
 * @author: T. Limpanuparb, J. Milthorpe
 */
public struct ShellPair {
    public val aang : Int; // Angular momentum at a
    public val bang : Int; // Angular momentum at b
    public val aPoint : Point3d; // Centre of a
    public val bPoint : Point3d; // Centre of b
    public val zetaA : Rail[Double]; // Gaussian exponent(s) on a
    public val zetaB : Rail[Double]; // Gaussian exponent(s) on b
    public val conA : Rail[Double]; // Contraction coefficient(s) for a
    public val conB : Rail[Double]; // Contraction coefficient(s) for b
    public val mu : Long;
    public val nu : Long;
    public val mu2 : Long;
    public val nu2 : Long;
    public val maxbraa : Int;
    public val maxbrab : Int;
    public val L : Rail[Int];
    public val maxL: Int; 
    public val contrib : Double;

    public def this(a:Int, b:Int, A:Point3d, B:Point3d, zetaA:Rail[Double], zetaB:Rail[Double], conA:Rail[Double], conB:Rail[Double], 
                    mu:Long, nu:Long, L:Rail[Int], contrib:Double) {
        this.aang=a;
        this.bang=b;
        this.aPoint=A;
        this.bPoint=B;
        this.zetaA=zetaA;
        this.zetaB=zetaB;
        this.conA=conA;
        this.conB=conB;
        this.mu=mu;
        this.nu=nu;
        this.maxbraa=(a+1n)*(a+2n)/2n;
        this.maxbrab=(b+1n)*(b+2n)/2n;
        this.mu2=mu+maxbraa-1n;
        this.nu2=nu+maxbrab-1n;
        this.L=L;        
        var maxl:Int=0n;
        for (i in (0..(L.size-1)))
           if (L(i)>maxl) maxl=L(i);
        this.maxL=maxl;
        this.contrib = contrib;
    }
}

