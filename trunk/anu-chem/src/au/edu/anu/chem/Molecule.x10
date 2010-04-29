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
    global val atomList = new ArrayList[T{self.at(this)}](); 
    global val name:String;

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

    public def getName() = this.name;

    public def addAtom(atm:T{self.at(this)}) : void {
        atomList.add(atm); 
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.i));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.j));
        maxExtent = Math.max(maxExtent, Math.abs(atm.centre.k));
    }

    public def getAtom(index:Int) : T{self.at(this)} = atomList.get(index) as T{self.at(this)};
    public def getAtoms() = atomList;
    public def getNumberOfAtoms() : Int = atomList.size();

    public def getNumberOfElectrons() : int {
       val ai = AtomInfo.getInstance();
       var ne:Int = 0;

       for(atm:T{self.at(this)} in atomList)
          ne += ai.getAtomicNumber(atm);

       return ne;
    }

    public def getMaxExtent() = maxExtent;

    public global safe def toString() : String {
       var str:String = "";

       for(atm:T{self.at(this)} in atomList)
         str += atm.toString() + "\n";

       return str;
    }
}

