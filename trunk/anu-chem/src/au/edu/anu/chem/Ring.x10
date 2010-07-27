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
package au.edu.anu.chem;

import x10.util.Pair;
import x10.util.ArrayList;
import x10.util.ValRailBuilder;

import x10x.vector.Point3d;

/**
 * This class represents a Ring structure in a Molecule
 *
 * @author: V.Ganesh
 */
public class Ring[T]{T <: Atom} {
    global val atomList = new ArrayList[T{self.at(this)}](); 
    var planar:Boolean;

    public def this() { 
        planar = true;   
    }

    public def addAtom(atm:T{self.at(this)}) : void {
        atomList.add(atm); 
    }

    public safe def getAtom(index:Int) : T{self.at(this)} = atomList.get(index) as T{self.at(this)};
    public global safe def getAtoms() = atomList;
    public global safe def getSize() = atomList.size();

    public def isPlanar() = planar;
    public def setPlanar(p:Boolean) { planar = p; }

    public global safe def toString() : String {
       var str:String = "[ ";

       for(atm:T{self.at(this)} in atomList)
         str += atm.toString() + " ";
       
       str += "]\n";

       return str;
    }
}

