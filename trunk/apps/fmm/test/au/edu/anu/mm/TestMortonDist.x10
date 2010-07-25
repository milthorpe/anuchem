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

import harness.x10Test;


/**
 * Test Morton distribution including generation and mapping to points of Morton indices.
 * @author milthorpe
 */
class TestMortonDist extends x10Test {
    public def run(): boolean {
        var boxRegion : Region(3) = [0..7, 0..7, 0..7];
        val boxDistribution = MortonDist.make(boxRegion);
        Console.OUT.println("boxDistribution: " + boxDistribution);

        val p = Point.make(0, 2, 1);
        val pm = MortonDist.getMortonIndex(p, 64);
        Console.OUT.println(pm.toBinaryString());
        Console.OUT.println(MortonDist.getPoint(pm, 64));

        Console.OUT.println(boxDistribution(p));

        val q = Point.make(1, 1, 0);
        val qm = MortonDist.getMortonIndex(q, 64);
        Console.OUT.println(qm.toBinaryString());
        Console.OUT.println(MortonDist.getPoint(qm, 64));
        Console.OUT.println(boxDistribution(q));

        val r = Point.make(5, 3, 7);
        val rm = MortonDist.getMortonIndex(r, 512);
        Console.OUT.println(rm.toBinaryString());
        Console.OUT.println(MortonDist.getPoint(rm, 512));
        Console.OUT.println(boxDistribution(r));

        return true;
    }

    public static def main(Rail[String]) {
        new TestMortonDist().execute();
    }

}
