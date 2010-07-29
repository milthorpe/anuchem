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

import x10.util.*;

/**
 * Orbitals.x10
 *
 * Represents an Orbital (used for basis function formation)
 *
 * @author: V.Ganesh
 */
public class Orbital { 
    global val exps:ArrayList[Double]{self.at(this)};
    global val coeff:ArrayList[Double]{self.at(this)};
    global val type:String{self.at(this)};
    global val angularMomentum:Int;

    public def this(type:String) { 
       exps      = new ArrayList[Double]();
       coeff     = new ArrayList[Double]();
       this.type = type as String{self.at(this)};

       if (type.equals("S")) angularMomentum = 0;
       else if (type.equals("P")) angularMomentum = 1;
       else if (type.equals("D")) angularMomentum = 2;
       else if (type.equals("F")) angularMomentum = 3;
       else angularMomentum = 0;
    } 

    public def addExponent(ex:Double) {
       exps.add(ex);
    }

    public def addCoefficient(co:Double) {
       coeff.add(co);
    }

    public def add(ex:Double, co:Double) {
       exps.add(ex);
       coeff.add(co);
    }

    public def getType() : String{self.at(this)} = this.type;
    public def getAngularMomentum() = this.angularMomentum;
    public def getExponents() : ArrayList[Double]{self.at(this)} = this.exps;
    public def getCoefficients() : ArrayList[Double]{self.at(this)} = this.coeff;
}

