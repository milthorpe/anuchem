package au.edu.anu.mm;

import au.edu.anu.chem.mm.MMAtom;

/**
 * This class represents a Morse potential between two atoms
 * of the form V(r) = D(1 - e^(-a(r-b)))^2 where:
 * D is the dissociation energy;
 * a is a constant related to the steepness of the curve (here a = 2/b); and
 * b is the equilibrium bond length
 * "The Morse curve is only a convenient analytical expression that has some
 * essential features of a diatomic potential, ... but there is no theoretical
 * justification for this particular form." - Berendsen, "Simulating the Physical World", p.6 (2007)
 * @see P. M. Morse "Diatomic molecules according to the wave mechanics.", Phys. Rev. 34, 57-64 (1929)
 */
public class DiatomicHarmonicPotential extends DiatomicPotential {
    /** The force constant in kJ mol^-1 nm^-2. */
    public val forceConstant : Double;

    /** The equilibrium bond length in nm. */
    public val bondLength : Double;

    public def this(atom1 : MMAtom, atom2 : MMAtom, bondLength : Double, forceConstant : Double) {
        super(atom1, atom2);
        this.bondLength = bondLength;
        this.forceConstant = forceConstant;
    }

    public def getPotentialAndForces() : Double {
        val r = atom2.centre - atom1.centre;
        val displacement = r.length() - bondLength;
        Console.OUT.println(r.length());
        atom1.force = forceConstant * displacement * r.normalize();
        atom2.force = -atom1.force;
        return 0.5 * forceConstant * displacement * displacement;
    }
}

