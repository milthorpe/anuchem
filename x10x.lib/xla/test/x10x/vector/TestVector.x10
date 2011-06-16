package x10x.vector;

import harness.x10Test;

/**
 * Tests XLA Vector class
 * @author milthorpe
 */
public class TestVector extends x10Test {
    public def run(): boolean {
        val vector = new Vector(3);
        for([i] in vector.vec) {
            vector.vec(i) = i as Double;
        }
        val magnitude : Double = vector.magnitude();
        chk(magnitude == Math.sqrt(5.0));
        return true;
    }

    public static def main(Array[String](1)) {
        new TestVector().execute();
    }
}
