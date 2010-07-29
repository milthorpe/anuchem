/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Australian National University 2010.
 */

package au.anu.edu.qm;

/**
 * Power.x10
 *
 * Represents a gaussian power
 *
 * @author: V.Ganesh
 */
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

