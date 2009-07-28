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
   
    var atomList:ArrayList[Atom];
    var name:String;

    public def this() { 
        atomList = new ArrayList[Atom]();
        name = "unknown";
    } 

    public def this(name:String) { 
        atomList  = new ArrayList[Atom]();
        this.name = name;
    } 

    public def getName() : String = this.name;
    public def addAtom(atm:Atom)  : void { atomList.add(atm); }
    public def getAtom(index:Int) : Atom = atomList.get(index);
    public def getAtoms() : ArrayList[Atom] = atomList;
    public def getNumberOfAtoms() : Int = atomList.size();

    public def setName(nam:String) : void { this.name = nam; }
    
    public def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public def toString() : String {
       var str:String = "";

       for(atm in atomList) 
         str += atm.getSymbol() + " " + atm.getX() + " " + atm.getY() + " " + atm.getZ() + "\n";

       return str;
    }
}

