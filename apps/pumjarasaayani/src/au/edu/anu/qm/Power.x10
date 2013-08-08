/*
 * This file is part of ANUChem.
 *
 *  This file is licensed to You under the Eclipse Public License (EPL);
 *  You may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *      http://www.opensource.org/licenses/eclipse-1.0.php
 *
 * (C) Copyright Australian National University 2010.
 */
package au.edu.anu.qm;

/**
 * Represents a gaussian power
 *
 * @author: V.Ganesh
 */
public struct Power(l:Int, m:Int, n:Int) { 
    val maxam:Int, minam:Int, totam:Int;

    public def this() { property(0n,0n,0n); maxam=minam=totam=0n; }

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

