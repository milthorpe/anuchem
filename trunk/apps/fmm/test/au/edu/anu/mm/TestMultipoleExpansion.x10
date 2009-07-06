package au.edu.anu.mm;

import harness.x10Test;
import x10x.vector.Point3d;

/**
 * Test multipole expansions
 * @author milthorpe
 */
class TestMultipoleExpansion extends x10Test {
    public def run(): boolean {
        val p : Int = 3; // multipole level

        val x : Point3d = new Point3d(1.0, 2.0, -1.0);
        val Olm : MultipoleExpansion = MultipoleExpansion.getOlm(1.5, x, p);
        Console.OUT.println("multipole expansion");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(Olm.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        
        val parent : Array[Complex]{rank==1} = Array.make[Complex]([0..(p+1)*(p+1)]->here, (val (i): Point)=> Complex.ZERO);
        var inm : Int = 0;
        for (val(n) : Point in [0..p]) {
            for (val(m) : Point in [-n..n]) {
                inm=inm+1;
                parent(inm) = Olm.terms(n,m);
            }
        }
        val child : Array[Complex]{rank==1} = Array.make[Complex]([0..(p+1)*(p+1)]->here, (val (i): Point)=> Complex.ZERO);

        MultipoleExpansion.translateExpansion(new Point3d(2.0, -3.0, 1.0), parent, child);
        Console.OUT.println("translated expansion");
		for (val(i) : Point in [0..(p+1)*(p+1)]) {
		    Console.OUT.print(child(i) + " ");
		}
        Console.OUT.println();

        val target : MultipoleExpansion = new MultipoleExpansion(p);
        MultipoleExpansion.translateAndAddMultipole(new Point3d(2.0, -3.0, 1.0), Olm, target);
        Console.OUT.println("translated multipole");
		for (val(i) : Point in [0..p]) {
            for (val(j) : Point in [-p..p]) {
			    Console.OUT.print(target.terms(i,j) + " ");
            }
            Console.OUT.println();
		}

        return true;
    }

    public static def main(Rail[String]) {
        new TestMultipoleExpansion().execute();
    }

}
