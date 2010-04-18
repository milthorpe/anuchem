package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

public abstract class DiatomicPotential {
    public val atom1 : MMAtom;
    public val atom2 : MMAtom;

    public def this(atom1 : MMAtom, atom2 : MMAtom) {
        this.atom1 = atom1;
        this.atom2 = atom2;
    }

    public abstract def getPotentialAndForces() : Double;
}

