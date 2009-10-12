/**
 * Molecule.x10
 *
 * Representation for a Molecule class 
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class Molecule { 
   
    global var atomList:ArrayList[Atom{self.at(this)}]{self.at(this)};
    global var name:String;

    public def this() { }

    public global def make() { 
        atomList = new ArrayList[Atom{self.at(this)}]();
        name = "unknown";
    } 

    public def make(name:String) { 
        atomList  = new ArrayList[Atom{self.at(this)}]();
        this.name = name;
    } 

    public def getName() : String = this.name;
    public def addAtom(atm:Atom{self.at(this)})  : void { atomList.add(atm); }
    public def getAtom(index:Int) : Atom{self.at(this)} = atomList.get(index) as Atom{self.at(this)};
    public def getAtoms() : ArrayList[Atom{self.at(this)}]{self.at(this)} = atomList;
    public def getNumberOfAtoms() : Int = atomList.size();

    public def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:Atom{self.at(this)} in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public def toString() : String {
       var str:String = "";

       for(atm:Atom{self.at(this)} in atomList) 
         str += atm.getSymbol() + " " + atm.getX() + " " + atm.getY() + " " + atm.getZ() + "\n";

       return str;
    }
}

