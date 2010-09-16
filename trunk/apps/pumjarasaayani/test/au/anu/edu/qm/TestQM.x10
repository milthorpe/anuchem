/**
 * Test QM code with various inputs
 * 
 * @author V. Ganesh
 */

package au.edu.anu.qm;

public class TestQM {
    public def this() {
    }

    private def runInput(inp:String, int:gMatMethod) : Double {
        val pm = new PumjaRasaayani(inp, gMatMethod);

        pm.runIt();

        return pm.getEnergy();
    }

    private def equalsTo(a:Double, b:Double, p:Double) : Boolean {
        return (Math.abs(a-b) < p);
    }

    public def run() {
        // iterate only once
        val jd = JobDefaults.getInstance();
        at(jd) { jd.setMaxIterations(1) };

        val path = "tests/";

        Console.OUT.println("Test result: " + equalsTo(runInput(path + "h2.inp", 1), -1.010695544230642, 1e-5));
    }

    public static def main(args : Array[String](1)) {
        new TestQM().run();
    } 
}

