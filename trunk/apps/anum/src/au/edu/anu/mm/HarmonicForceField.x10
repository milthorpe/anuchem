package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class HarmonicForceField implements ForceField {
    val diatomicPotentials : ValRail[DiatomicHarmonicPotential];

    public def this(diatomicPotentials : ValRail[DiatomicHarmonicPotential]) {
        this.diatomicPotentials = diatomicPotentials;
    }
    
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double {
        var V : Double = 0.0;        
        for((p) in 0..diatomicPotentials.length-1) {
            val potential = diatomicPotentials(p);
            V += potential.getPotentialAndForces();
        }
        return V;
    }
}
