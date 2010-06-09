/**
 * Ring.x10
 * This class represents a Ring structure in a Molecule
 *
 * @author: V.Ganesh
 */
package au.edu.anu.chem;

import x10.util.Pair;
import x10.util.ArrayList;
import x10.util.ValRailBuilder;

import x10x.vector.Point3d;

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

