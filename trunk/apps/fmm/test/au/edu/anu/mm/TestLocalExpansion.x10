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
import x10x.vector.Tuple3d;

/**
 * Test local Taylor-type expansions.
 * @author milthorpe
 */
class TestLocalExpansion extends MathTest {
    public def run(): boolean {
        val p = 3; // multipole level
        val x = new Point3d(1.0, 2.0, -1.0);

        val Mlm = LocalExpansion.getMlm(x as Tuple3d, p);
        Console.OUT.println("local expansion:\n" + Mlm.toString());

        val target = new LocalExpansion(p);
        val y = new Point3d(2.0, -3.0, 1.0);
        val translation = MultipoleExpansion.getOlm(y as Tuple3d, p);
        target.translateAndAddLocal(translation, Mlm);
        Console.OUT.println("translated local:\n" + target.toString());

        val roundtrip = new LocalExpansion(p);
        val z = new Point3d(-2.0, 3.0, -1.0);
        val reverseTranslation = MultipoleExpansion.getOlm(z as Tuple3d, p);
        roundtrip.translateAndAddLocal(reverseTranslation, target);
        Console.OUT.println("translated local - roundtrip:\n" + roundtrip.toString());
		for ((i) in 0..p) {
            for ((j) in -i..i) {
                chk(nearlyEqual(roundtrip.terms(i,j), Mlm.terms(i,j), 1.0e-6, 1.0e-12));
		    }
        }

        return true;
    }

    public static def main(Rail[String]) {
        new TestLocalExpansion().execute();
    }

}
