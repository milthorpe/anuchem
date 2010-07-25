/*
 * This file is part of ANUChem.
 * ANUChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * ANUChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with ANUChem.  If not, see <http://www.gnu.org/licenses/>.
 *
 * (C) Copyright Josh Milthorpe 2010.
 */
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

