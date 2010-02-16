/**
 * This class represents an Atom for the purpose of quantum
 * chemical and other computational chemistry codes.
 *
 * @author milthorpe
 *
 * @modified ganesh - toString() to simply print coordinate values insted of 
 *                    default (i, j , k) print.
 */
package au.edu.anu.chem;

import x10x.vector.Point3d;

public class Atom { 
    /** The location of the atomic nucleus. */
    public global var centre : Point3d!;

    /** The symbol for this atom species. */
    public global val symbol : String;

    public def this(symbol : String, centre : Point3d!) { 
        this.symbol = symbol;
        this.centre = centre;
    }

    public def this(centre : Point3d!) { 
        this.symbol = "";
        this.centre = centre;
    }

    public global safe def toString() : String {
        return symbol + " " + centre.i + " " + centre.j + " " + centre.k;
    }
}

