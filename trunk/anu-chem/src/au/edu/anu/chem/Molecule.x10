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

    public def getName() : String = this.name;
    public def addAtom(atm:T{self.at(this)})  : void { atomList.add(atm); }
    public def getAtom(index:Int) : T{self.at(this)} = atomList.get(index) as T{self.at(this)};
    public def getAtoms() : ArrayList[T{self.at(this)}]{self.at(this)} = atomList;
    public def getNumberOfAtoms() : Int = atomList.size();

    public def getNumberOfElectrons() : int {
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

