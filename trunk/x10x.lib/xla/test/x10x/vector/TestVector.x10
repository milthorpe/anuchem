package x10x.vector;

import harness.x10Test;

/**
 * Tests XLA Vector class
 * @author milthorpe
 */
public class TestVector extends x10Test {
    public def run(): boolean {
        val vector = Vector.make(3) as Vector!;
        ateach(p(i) : Point in vector.vec) {
            vector.vec(i) = i as Double;
        }
        // TODO magnitude not place-safe
        //val magnitude : Double = vector.magnitude();
        // TODO Vector.make returns vector of wrong size!
        //chk(magnitude == Math.sqrt(5.0));
        return true;
    }

    public static def main(Array[String](1)) {
        new TestVector().execute();
    }
}
