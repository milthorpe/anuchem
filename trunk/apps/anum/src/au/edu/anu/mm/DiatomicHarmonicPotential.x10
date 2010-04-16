package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

public class DiatomicHarmonicPotential {
    public val atom1 : MMAtom;
    public val atom2 : MMAtom;

    /** The force constant in kJ mol^-1 nm^-2. */
    public val forceConstant : Double;

    /** The equilibrium bond length in nm. */
    public val bondLength : Double;

    public def this(atom1 : MMAtom, atom2 : MMAtom, bondLength : Double, forceConstant : Double) {
        this.atom1 = atom1;
        this.atom2 = atom2;
        this.bondLength = bondLength;
        this.forceConstant = forceConstant;
    }

    public def getPotentialAndForces() : Double {
        val r = atom2.centre - atom1.centre;
        val displacement = r.length() - bondLength;
        atom1.force = forceConstant * displacement * r.normalize();
        atom2.force = -atom1.force;
        return 0.5 * forceConstant * displacement * displacement;
    }
}

