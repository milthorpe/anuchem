package au.edu.anu.mm;

import x10x.vector.Point3d;
import au.edu.anu.chem.mm.MMAtom;
import au.edu.anu.chem.mm.ElectrostaticDirectMethod;
import au.edu.anu.chem.mm.TestElectrostatic;
import au.edu.anu.util.Timer;

/**
 * Tests a simple pairwise electrostatic calculation.
 * @author milthorpe
 */
public class TestPairwise extends TestElectrostatic {
    public global def sizeOfCentralCluster() : Double = 80.0;

    public def this(numAtoms : Int) {
        super(numAtoms);
    }

    public static def main(args : Rail[String]!) {
        var numAtoms : Int;
        if (args.length > 0) {
            numAtoms = Int.parseInt(args(0));
        } else {
            Console.ERR.println("usage: TestFmm3d numAtoms [density] [numTerms] [wellSpaced]");
            return;
        }

        new TestPairwise(numAtoms).test();
    }

    public def test() {
        Console.OUT.println("Testing pairwise electrostatic calculation for " + numAtoms + " atoms");

        val atoms = generateAtoms();
        val direct = new ElectrostaticDirectMethod(atoms);
        val directEnergy = direct.getEnergy();
        logTime("Total time", ElectrostaticDirectMethod.TIMER_INDEX_TOTAL, direct.timer);
    }
}

