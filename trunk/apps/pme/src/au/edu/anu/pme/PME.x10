package au.edu.anu.pme;

import x10x.vector.Point3d;

public class PME {
	/** The cartesian location of the top-left-front corner of the simulation cube. */
    public val topLeftFront : Point3d;

    /** The length of a side of the simulation cube. */
    public val size : Double; 

	val atoms : ValRail[Atom];
	
	public def this(topLeftFront : Point3d,
            size : Double,
            atoms : ValRail[Atom]) {
		this.topLeftFront = topLeftFront;
        this.size = size;

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
	
	
}