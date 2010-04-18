package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

public interface ForceField {
    public global def getPotentialAndForces(atoms: DistArray[ValRail[MMAtom]](1)) : Double;
}
