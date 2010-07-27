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
 * This class represents a Molecule
 * for the purposes of computational chemistry codes.
 *
 * @author: V.Ganesh
 */
public class Molecule[T]{T <: Atom} {
    global val atomList = new ArrayList[T{self.at(this)}](); 
    global val name : String;

    global val ringList = new ArrayList[Ring[T]{self.at(this)}]();

    /** 
     * Measures the maximum absolute value of any coordinate x,y,z
     * of all atoms. This is used to estimate a rough cubic box size.
     */
    private var maxExtent : Double = 0.0;

    public def this() { 
        name = "unknown";
    }

    public def this(name:String) {
        this.name = name;
    }

    public global safe def getName() = this.name;

    public def addAtom(atm:T{self.at(this)}) : void {
        atomList.add(atm); 
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.i));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.j));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.k));      
    }

    public safe def getAtom(index:Int) : T{self.at(this)} = atomList.get(index) as T{self.at(this)};
    public global safe def getAtoms() = atomList;
    public global safe def getNumberOfAtoms() : Int = atomList.size();

    public global safe def addRing(r:Ring[T]{self.at(this)}) { ringList.add(r); }
    public global safe def getRings() = ringList;

    public safe def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:T{self.at(this)} in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public safe def getMaxExtent() = maxExtent;

    public safe def getCoords() : ValRail[Pair[String,Point3d]] {
        val coords = new ValRailBuilder[Pair[String,Point3d]](atomList.size());
        for(atom in atomList) {
            coords.add(Pair[String,Point3d](atom.symbol, atom.centre));
        }
        return coords.result();
    }

    public safe def centreOfMass() : Point3d {
        val ai = AtomInfo.getInstance();
        var x:Double = 0.0, y:Double = 0.0, z:Double = 0.0;
        var massSum:Double = 0.0;

        for(atm:T{self.at(this)} in atomList) {
            val mass = ai.getAtomicMass(atm);
            x += mass * atm.centre.i;
            y += mass * atm.centre.j;
            z += mass * atm.centre.k;

            massSum += mass;
        }

        return Point3d(x/massSum,y/massSum,z/massSum);
    }

    public global safe def toString() : String {
       var str:String = "";

       for(atm:T{self.at(this)} in atomList)
         str += atm.toString() + "\n";

       return str;
    }
}

