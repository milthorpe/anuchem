/**
 * Atom.x10
 *
 * Represent an Atom class
 *
 * @author: V.Ganesh
 */

package au.anu.edu.qm;

public value Atom { 
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

    public def getX() : Double = this.x; 
    public def getY() : Double = this.y; 
    public def getZ() : Double = this.z; 

    public def getSymbol() : String = this.symbol; 

    public def distanceFrom(atm:Atom) : Double {
        return Math.sqrt(distanceSquaredFrom(atm));
    }

    public def distanceSquaredFrom(atm:Atom) : Double {
        val x = this.x - atm.x;
        val y = this.y - atm.y;
        val z = this.z - atm.z;

        return (x*x + y*y + z*z);
    }
}

