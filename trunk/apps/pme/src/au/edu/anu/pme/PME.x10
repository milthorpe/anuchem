package au.edu.anu.pme;

import x10x.vector.Point3d;
import x10x.vector.Vector3d;

public class PME {
	/** The cartesian location of the top-left-front corner of the simulation cube. */
    public val topLeftFront : Point3d;

    /** The length of a side of the simulation cube. */
    public val size : Double;

    /** The edges of the unit cell. TODO: non-cubic cells */
    public val edges : ValRail[Vector3d](3);

    /** The conjugate reciprocal vectors for each dimension. */
    public val edgeReciprocals : ValRail[Vector3d](3);

	val atoms : ValRail[Atom];
	
    public def this(topLeftFront : Point3d,
            size : Double,
            atoms : ValRail[Atom]) {
	    this.topLeftFront = topLeftFront;
        this.size = size;
        this.edges = [new Vector3d(size, 0.0, 0.0), new Vector3d(0.0, size, 0.0), new Vector3d(0.0, 0.0, size)];
        this.edgeReciprocals = [new Vector3d(1.0 / size, 0.0, 0.0), new Vector3d(0.0, 1.0 / size, 0.0), new Vector3d(0.0, 0.0, 1.0 / size)];
        this.atoms = atoms;
    }
	
    public def calculateEnergy() : Double {
        var directEnergy : Double = 0.0;
        for ((i) in 0..(atoms.length - 1)) {
            for ((j) in 0..(i - 1)) {
                val pairEnergy : Double = atoms(j).pairEnergy(atoms(i));
                directEnergy += 2 * pairEnergy;
            }
        }
         Console.OUT.println("directEnergy = " + directEnergy);
         return directEnergy;
    }
    
    public def getScaledFractionalCoordinates(r : Vector3d) : Vector3d {
        return new Vector3d(r.dot(edgeReciprocals(0)), r.dot(edgeReciprocals(1)), r.dot(edgeReciprocals(2)));
    }


    public def w2(u : Double) = Math.abs(u) <= 1 ? 1 : 0;
}
