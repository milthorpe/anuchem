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