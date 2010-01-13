/**
 * Atom.x10
 *
 * Represent an Atom class
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

import x10.util.*;

public class Atom { 
    global val x:Double, y:Double, z:Double;
    global val symbol:String;

    public def this(symbol:String, x:Double, y:Double, z:Double) { 
        this.symbol = symbol;
        this.x = x; this.y = y; this.z = z;
    } 

    public def this(x:Double, y:Double, z:Double) { 
        this.symbol = "";
        this.x = x; this.y = y; this.z = z;
    } 

    public def getX() = this.x; 
    public def getY() = this.y; 
    public def getZ() = this.z; 

    public def getSymbol() = this.symbol; 

    
    global var basisFunctions:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)};
    
    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)}) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() : ArrayList[ContractedGaussian{self.at(this)}]{self.at(this)} = basisFunctions;    

    public def distanceFrom(atm:Atom) : Double {
        return Math.sqrt(distanceSquaredFrom(atm));
    }

    public def distanceSquaredFrom(atm:Atom) : Double {
        val x = this.x - atm.x;
        val y = this.y - atm.y;
        val z = this.z - atm.z;

        return (x*x + y*y + z*z);
    }

    public global safe def toString() : String {
        return symbol + " " + x + " " + y + " " + z;
    }
}

