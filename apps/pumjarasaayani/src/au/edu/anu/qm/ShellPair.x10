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
    public val dconA : Int; // Degree of contraction on a
    public val dconB : Int; // Degree of contraction on b
    public val mu : Long;
    public val nu : Long;
    public val mu2 : Long;
    public val nu2 : Long;
    public val maxbraa : Int;
    public val maxbrab : Int;
    public val maxL : Rail[Long];
    public val contrib : Double;

    public def this(a:Int, b:Int, A:Point3d, B:Point3d, zetaA:Rail[Double], zetaB:Rail[Double], conA:Rail[Double], conB:Rail[Double], 
                    dconA:Int, dconB:Int, mu:Long, nu:Long, maxL:Rail[Long], contrib:Double) {
        this.aang=a;
        this.bang=b;
        this.aPoint=A;
        this.bPoint=B;
        this.zetaA=zetaA;
        this.zetaB=zetaB;
        this.conA=conA;
        this.conB=conB;
        this.dconA=dconA; 
        this.dconB=dconB; 
        this.mu=mu;
        this.nu=nu;
        this.maxbraa=(a+1)*(a+2)/2;
        this.maxbrab=(b+1)*(b+2)/2;
        this.mu2=mu+maxbraa-1;
        this.nu2=nu+maxbrab-1;
        this.maxL=maxL;        
        this.contrib = contrib;
    }
}

