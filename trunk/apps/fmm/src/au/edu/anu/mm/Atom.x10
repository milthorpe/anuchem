package au.edu.anu.mm;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

/**
 * This class represents an Atom for the purpose of molecular 
 * mechanics calculations.
 * // TODO connectivity for force fields
 * @author milthorpe
 */
public class Atom { 
    /** The location of the atomic nucleus. */
    public var centre : Point3d;

    /** The current force acting upon this atom. */
    public var force : Vector3d;

    /** The atomic symbol e.g. "Li" */
    public val symbol : String;

    /** The effective charge in atomic units. */
    public val charge : Double;

    public def this(symbol : String, centre : Point3d, charge : Double) { 
        this.symbol = symbol;
        this.centre = centre;
        this.charge = charge;
        this.force = Vector3d.NULL;
    } 

    public def this(symbol:String, centre : Point3d) { 
        this(symbol, centre, 0.0);
    } 
}

