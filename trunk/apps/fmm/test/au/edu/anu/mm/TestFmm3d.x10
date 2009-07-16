package au.edu.anu.mm;

import x10x.vector.Point3d;

/**
 * Tests the FMM 3D implementation.
 * @author milthorpe
 */
public class TestFmm3d {
    public static def main(args: Rail[String]) {
        var numLevels : Int;
        if (args.length >= 1) {
            numLevels = Int.parseInt(args(0));
        } else {
            numLevels = 3;
            //Console.ERR.println("usage: TestFmm3d <numLevels>");
        }
        val atoms : Rail[Atom] = ValRail.make[Atom](1, (Int) => new Atom("H", new Point3d(0, 0, 0)));
        val energy : Double = new Fmm3d(numLevels, atoms).calculateEnergy();
    }
}

