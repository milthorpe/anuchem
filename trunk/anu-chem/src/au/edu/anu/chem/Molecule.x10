/**
 * Molecule.x10
 * This class represents a Molecule
 * for the purposes of computational chemistry codes.
 *
 * @author: V.Ganesh
 */
package au.edu.anu.chem;

import x10.util.ArrayList;

public class Molecule[T]{T <: Atom} {
    global var atomList:ArrayList[T{self.at(this)}]{self.at(this)};
    global var name:String;

    public def this() { }

    public global def make() {
        atomList = new ArrayList[T{self.at(this)}]();
        name = "unknown";
    }

    public def make(name:String) {
        atomList  = new ArrayList[T{self.at(this)}]();
        this.name = name;
    }

    public global def getName() : String = this.name;
    public global def addAtom(atm:T{self.at(this)})  : void { atomList.add(atm); }
    public global def getAtom(index:Int) : T{self.at(this)} = atomList.get(index) as T{self.at(this)};
    public global def getAtoms() : ArrayList[T{self.at(this)}]{self.at(this)} = atomList;
    public global def getNumberOfAtoms() : Int = atomList.size();

    public global def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:T{self.at(this)} in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public global safe def toString() : String {
       var str:String = "";

       for(atm:T{self.at(this)} in atomList)
         str += atm.toString() + "\n";

       return str;
    }
}

