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
    val x:Double, y:Double, z:Double;
    val symbol:String;

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

    
    var basisFunctions:ArrayList[ContractedGaussian];
    
    public def setBasisFunctions(bfs:ArrayList[ContractedGaussian]) {
        basisFunctions = bfs;
    }

    public def getBasisFunctions() = basisFunctions;    

    public def distanceFrom(atm:Atom) {
        return Math.sqrt(distanceSquaredFrom(atm));
    }

    public def distanceSquaredFrom(atm:Atom) {
        val x = this.x - atm.x;
        val y = this.y - atm.y;
        val z = this.z - atm.z;

        return (x*x + y*y + z*z);
    }
}

