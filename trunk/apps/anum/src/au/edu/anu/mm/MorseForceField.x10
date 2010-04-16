package au.edu.anu.mm;

import x10.util.Pair;
import x10x.vector.Vector3d;
import au.edu.anu.chem.mm.MMAtom;

public class MorseForceField implements ForceField {
    val diatomicMolecules : ValRail[Pair[MMAtom,MMAtom]];

    public def this(diatomicMolecules : ValRail[Pair[MMAtom,MMAtom]]) {
        this.diatomicMolecules = diatomicMolecules;
    }
    
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double {
        finish ateach((p) in atoms) {
            val myAtoms = atoms(p);
            finish foreach((i) in 0..myAtoms.length-1) {
                myAtoms(i).force = new Vector3d(0.1, 0.0, 0.0);
            }
        }
        return 0.0;
    }
}
